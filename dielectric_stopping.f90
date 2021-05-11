module DielectricStopping
  use utils, only : dp, pi
  use RPAdielectric, only : FDD, elf
  use MerminDielectric, only : mELF
  ! sorting package from MJ Rutter  
  use sorts, only         : quicksort
  use random, only        : normal
  use integrate, only     : trapezoidal
  use ArrayOps, only      : linspace, logspace, bunchspace, binarysearch
  
  implicit none
  private
  public omegaintegral, stoppingintegral, sumrule, srdensity, FEGdensity, wp, BG_wp, modBG_wp 
  
contains
  
  ! subroutine sortpositive(a, inp)
  !   ! Sort numbers in array a (in place) and note the index that separates the
  !   ! negative points from the positive points
  !   real(dp), intent(inout) :: a(:)
  !   integer, intent(out) :: inp
  !   integer :: len, i

  !   len = size(a)
  !   call quicksort(a)

  !   do i = 1,len
  !      if (a(i) >= 0.) then
  !         inp = i
  !         exit
  !      end if
  !   end do
    
  ! end subroutine sortpositive
    
  function sumrule(n)
    real(dp), intent(in) :: n
    real(dp) :: sumrule

    sumrule = pi/2 * wp(n)**2

    return 
  end function sumrule

  function FEGdensity(kbT, mu) result(n)
    ! Calculate density of a free-electron gas given the temperature and
    ! chemical potential. Everything in atomic units.
    ! Currently assuming a smooth integrand, i.e. relatively large temperatures.

    real(dp), intent(in) :: kbT, mu
    real(dp) :: n
    integer  :: i
    real(dp) :: x(201)
    real(dp) :: y(size(x))

    ! ! integrate
    x = linspace(0._dp, mu + 10*kbT, size(x))!bunchspace(mu, mu - mu*1e-4, (size(x)+1)/2)
    do i=1,size(x)
       y(i) = FDD(x(i), mu, kbT, sqrt(x(i)))
    end do
    n = trapezoidal(y,x)
    
    n = sqrt(2._dp)/pi**2 * n
    
  end function FEGdensity

  function srdensity(kbT, mu) result(den)
    !! Compute density using the f-sum rule
    real(dp), intent(in) :: kBT, mu
    real(dp) :: den
    real(dp) :: k, v, n, int ! dummy variables
    logical  :: srcheck

    k = 1.
    ! Need an approximation to the density to start
    n = FEGdensity(kbT, mu)

    srcheck = .false.
    do while (.not. srcheck)
       ! In omegaintegral, we want k*v > ELFmaxpos + ELFwidth
       v = modBG_wp(n, kbT, mu)/k + sqrt(2*(20*kbT + mu))
       call omegaintegral(v, k, kbT, mu, n, int, srcheck)
       k = k + 1.
    end do

    ! int should be the sum rule
    den = int / (2 * pi**2)
   
  end function srdensity
  
  function wp(n)
    ! Basic plasma frequency, a.u.
    real(dp), intent(in) :: n
    real(dp) :: wp

    wp = sqrt(4*pi*n)
  end function wp

  function BG_wp(n, k, kbT)
    ! Bohm-Gross dispersion relation for the plasma frequency, in a.u.

    real(dp), intent(in) :: n, k, kbT
    real(dp) :: BG_wp

    BG_wp = sqrt(wp(n)**2 + 3*kbT*k**2)
  end function BG_wp

  function modBG_wp(n, k, kbT)
    ! Modified Bohm-Gross dispersion relation, given in Glezner & Redmer, Rev.
    ! Mod. Phys., 2009, Eqn (16).
    real(dp), intent(in) :: n, k, kbT
    real(dp) :: modBG_wp, therm_deBroglie
    therm_deBroglie = sqrt(2* pi / kbT)

    modBG_wp = sqrt(BG_wp(n, k, kbT)**2 + 3*kbT*k**2 &
               + 0.088 * n * therm_deBroglie**3 + k**4/4)
  end function modBG_wp

  subroutine omegaintegral(v, k, kbT, mu, density, omegaint, SRsatisfy)
    real(dp), intent(in)    :: v, k, kbT, mu, density
    real(dp), intent(inout) :: omegaint
    logical , intent(inout) :: SRsatisfy
    real(dp) :: int_allspace ! integral over all space, should give sum rule
    real(dp) :: sr ! sum rule value
    real(dp) :: ELFwidth, ELFmaxpos, ELFmin, ELFmax ! ELF variables
    ! real(dp), allocatable :: w(:), y(:) ! integration variable and integrand
    ! Three ranges in omega-space corresponding to 3 possible relative
    ! locations of k*v (upper limit of omega interval)
    real(dp) :: wminrange(100), wmaxrange(100), wELFrange(999)
    real(dp) :: yminrange(size(wminrange)), ymaxrange(size(wmaxrange))
    real(dp) :: yELFrange(size(wELFrange))
    real(dp) :: ytemp, wtemp ! tempararily holds values from w, y arrays.
    integer  :: i, N ! Number of steps for integration

    ! sum rule value
    sr = sumrule(density)

    ! Max position of ELF peak (approximation)
    ELFmaxpos = modBG_wp(density, k, kbT)
    ! approximate width of ELF peak
    ELFwidth = sqrt(2*(20*kbT + mu))*k
    ! Leftmost edge of ELF peak
    ELFmin = max(0.001_dp, ELFmaxpos - ELFwidth)
    ! Want ELFwidth < ELFmaxpos
    ELFwidth = ELFmaxpos - ELFmin
    ! rightmost edge of ELF peak
    ELFmax = ELFmaxpos + ELFwidth

    if (k*v < ELFmin) then
       ! kv < ELFmin: ELF peak not in integration range

       wminrange = linspace(0.001_dp, k*v, size(wminrange))
       do i = 1, size(wminrange)
          yminrange(i) = wminrange(i) * mELF(k, wminrange(i), kbT, mu)!, &
               ! 1e-3_dp, 0._dp)
       end do
       
       omegaint = trapezoidal(yminrange, wminrange)
       ! Don't check in this case since we are not dealing
       ! with the difficult part (peak) of the integrand
       SRsatisfy = .true. 
       
       return
    end if
    
    ! ELFmin < kv < ELFmax or ELFmax < kv

    ! integral from (0, ELFmin)
    wminrange = linspace(0.001_dp, ELFmin, size(wminrange))
    do i = 1, size(wminrange)
       yminrange(i) = wminrange(i) * mELF(k, wminrange(i), kbT, mu!), &
            !1e-3_dp, 0._dp)
    end do
    
    omegaint = trapezoidal(yminrange, wminrange)
    
    ! points within the ELF peak
    wELFrange = bunchspace(ELFmaxpos, ELFwidth, (size(wELFrange)+1)/2)
    
    ! compute integrand for each point in w
    do i = 1, size(wELFrange)
       yELFrange(i) = wELFrange(i) * mELF(k, wELFrange(i), kbT, mu!), &
            !1e-3_dp, 0._dp)
    end do
    
    ! ELFmin < kv < ELFmax
    if ((ELFmin < k*v) .and. (k*v < ELFmax)) then
       ! ! integral from (ELFmin, kv) and (kv, ELFmax)

       ! Find index corresponding to the upper limit of the integral (k*v)
       i = binarysearch(wELFrange, k*v)
       ! integral from [ELFmin, k*v]
       wtemp = wELFrange(i)
       ytemp = yELFrange(i)
       wELFrange(i) = k*v
       yELFrange(i) = wELFrange(i) * mELF(k, wELFrange(i), kbT, mu!), &
            !1e-3_dp, 0._dp)
       omegaint = omegaint + trapezoidal(yELFrange(1:i), wELFrange(1:i))
       
       ! integral from [kv, ELFmax \approx \infty]
       ! (for integration over (0, \infty))
       wELFrange(i-1) = wELFrange(i)
       wELFrange(i) = wtemp
       yELFrange(i-1) = yELFrange(i)
       yELFrange(i) = ytemp
       int_allspace = omegaint + &
            trapezoidal(yELFrange(i-1:size(yELFrange)), &
                        wELFrange(i-1:size(wELFrange)))
       
    else ! kv > ELFmax
       ! integral from (ELFmin, ELFmax)
       omegaint = omegaint + trapezoidal(yELFrange, wELFrange)
       ! integral from (ELFmax, kv)
       wmaxrange = linspace(ELFmax, k*v, size(wmaxrange))
       do i = 1, size(wmaxrange)
          ymaxrange(i) = wmaxrange(i) * mELF(k, wmaxrange(i), kbT, mu!), &
               !1e-3_dp, 0._dp)
       end do
       
       omegaint = omegaint + trapezoidal(ymaxrange, wmaxrange)
       ! In this case, the integral over (0, \infty) is the same as the
       ! integral over (0, kv)
       int_allspace = omegaint
    end if

    ! Is the sum rule satisfied?
    if (abs(sr - int_allspace) / sr > 1e-2)  then
       ! print *, "Convergence to sum rule failed in omega integral:"
       ! print "(a, e7.2)", "k = ", k
       ! print "(a, f9.4)", "sum rule = ", sr
       ! print "(a, f9.4)", "integral = ", int_allspace
       ! print "(a, f9.4)", "integral up to kv = ", omegaint
       SRsatisfy = .false.
    else
       SRsatisfy = .true.
    end if

  end subroutine omegaintegral

  subroutine stoppingintegral(v, kbT, mu, kint, upper_bound_err)
    real(dp), intent(in)    :: v, kbT, mu
    real(dp), intent(inout) :: kint ! stopping power integral
    real(dp), intent(inout) :: upper_bound_err ! upper bound on error
    real(dp)                :: density, width, kint_upperbound
    real(dp)                :: k(200) ! integration variable
    real(dp)                :: y(size(k)) ! integrand
    integer                 :: i, N ! Number of steps for integration
    logical                 :: SRsatisfy
    real(dp)                :: klowest ! Lowest value of k for which we satisfy
                                       ! the sum rule
    
    density = srdensity(kbT, mu)

    N = size(k)
    width = (2*(10*kbT + mu))**(0.5)
    kint_upperbound = 2*(v +  width)
    k = logspace(0.01_dp, kint_upperbound, N)
    y = 0. ! sets all elements to 0
    upper_bound_err = 0.
    klowest = 0.
    SRsatisfy = .true.
    ! call omegaintegral(v, k(5), kbT, mu, density, y(5), SRsatisfy)
    !print *,  y(1)
    do i = 1, N !N, 1, -1
       ! print *, i
       call omegaintegral(v, k(i), kbT, mu, density, y(i), SRsatisfy)
       ! print *, SRsatisfy
       y(i) = y(i) / k(i)
       if (.not. SRsatisfy) then
          ! If sum rule is not satisfied for some k*, it is unlikely to be
          ! satisfied for any k <= k* (the omega integrand becomes
          ! progressively harder to integrate).
          ! print *, "SR not satisfied, v = ", v, " k = ", k(i)
          klowest = k(i)
          y(i) = 0
          ! exit
       end if
    end do
    
    kint = trapezoidal(y, k)
    if (klowest * v > wp(density)) upper_bound_err = pi/2 * klowest * v
    
  end subroutine stoppingintegral
  
end module DielectricStopping
