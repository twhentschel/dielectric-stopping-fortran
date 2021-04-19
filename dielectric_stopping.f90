module DielectricStopping
  use RPAdielectric, only : dp, pi, FDD, elf
  ! sorting package from MJ Rutter  
  use sorts, only         : quicksort
  use random, only        : normal
  use integrate, only     : trapezoidal
  use fortranArrays, only : linspace, bunchspace
  
  implicit none
  
  public
  
contains
  
  subroutine sortpositive(a, inp)
    ! Sort numbers in array a (in place) and note the index that separates the
    ! negative points from the positive points
    real(dp), intent(inout) :: a(:)
    integer, intent(out) :: inp
    integer :: len, i

    len = size(a)
    call quicksort(a)

    do i = 1,len
       if (a(i) >= 0.) then
          inp = i
          exit
       end if
    end do
    
  end subroutine sortpositive
    
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
    real(dp) :: n, x, xmax, dx
    integer  :: i, imax

    xmax = mu + 10*kbT ! 1/exp(10) \approx 1e-9
    imax = 200
    dx = xmax/imax

    ! integrate
    n = 0.0
    do i=1,imax
       x = i*dx
       n = n + FDD(x, mu, kbT, sqrt(x))*dx
    end do

    n = sqrt(2._dp)/pi**2 * n
    
  end function FEGdensity

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
  
  function omegaint(v, k, kbT, mu, ne)
    real(dp), intent(in) :: v, k, kbT, mu
    real(dp), intent(in), optional :: ne ! density
    real(dp) :: omegaint, omegaint2, intmin, intmax, dw, density, sr
    real(dp) :: ELFwidth, ELFwmax, ELFmin, ELFmax
    real(dp), allocatable :: w(:), y(:) ! integration variable and integrand
    integer  :: i, ipos ! do-while loop
    integer  :: N ! Number of steps for integration

    if (present(ne)) then
       density = ne
    else
       density = FEGdensity(kbT, mu)
    end if

    ! sum rule value
    sr = sumrule(density)

    ! Max position of ELF peak (approximation)
    ELFwmax = modBG_wp(density, k, kbT)
    ! approximate width of ELF peak
    ELFwidth = sqrt(2*(20*kbT + mu))*k
    ! Leftmost edge of ELF peak
    ELFmin = max(0.001_dp, ELFwmax - ELFwidth)
    ELFwidth = ELFwmax - ELFmin
    ! rightmost edge of ELF peak
    ELFmax = ELFwmax + ELFwidth

    print *, "wp = ", wp(density)
    print *, "kv = ", k*v
    print *, "ELFwidth = ", ELFwidth
    print *, "ELfmin = ", ELFmin
    print *, "Elfwmax = ", ELFwmax
    print *, "ELFmax = ", ELFmax

    ! N = 100
    ! allocate(w(N))

    ! call normal(w)
    ! w = ELFwidth * w + ELFwmax
    ! call sortpositive(w, ipos)
    ! allocate(y(N-(ipos-1)))

    N = 1000
    allocate(w(2*N-1))
    allocate(y(2*N-1))
    w = bunchspace(ELFwmax, ELFwidth, N)

    ! compute integrand for each point in w
    do i = 1, size(w)
       y(i) = w(i) * elf(k, w(i), kbT, mu)
    end do

    omegaint = trapezoidal(y, w)

    ! ! compute integrand for each point in w
    ! do i = ipos, N
    !    y(i-(ipos-1)) = w(i) * elf(k, w(i), kbT, mu)
    ! end do

    ! omegaint = trapezoidal(y, w(ipos:N))
    
    deallocate(y) 
    deallocate(w)
    
    
    print "(a, f9.6)", "sumrule = ", sr
    print "(a, f9.6)", "integral = ", omegaint
    
    ! ! Nothing in integration range
    ! if (k*v < ELFmin) then
    !    omegaint = 0.
    !    return
    ! end if

    ! ! integrate first part
    ! intmax = min(k*v, ELFmax)
    ! intmin = ELFmin
    ! imax = 1500
    ! dw = (intmax - intmin)/imax

    ! omegaint = 0.
    ! do i=1,imax
    !    w = i*dw + intmin
    !    omegaint = omegaint + w * elf(k, w, kbT, mu) * dw
    ! end do

    ! ! integrate the second part
    ! if (intmax == ELFmax) then
    !    omegaint2 = 0.
    ! else
    !    intmax = ELFmax
    !    intmin = k*v
    !    imax = 500
    !    dw = (intmax - intmin)/imax

    !    omegaint2 = 0.
    !    do i=1,imax
    !       w = i*dw + intmin
    !       omegaint2 = omegaint2 + w * elf(k, w, kbT, mu) * dw
    !    end do
    ! end if

    ! omegaint = omegaint + omegaint2

    ! print *, omegaint2
    
  end function omegaint

  function kint(v, kbT, mu, ne)
    
  end function kint
end module DielectricStopping
