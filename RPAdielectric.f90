module RPAdielectric
  use utils, only : dp, pi
  implicit none
  private
  public FDD, redielectric, imdielectric, elf
  
contains
  
  function FDD(E, a, b, g_in)
    ! Fermi-Dirac Distriution function.
    real(dp), intent(in) :: E, a, b
    real(dp), intent(in), optional :: g_in ! density of states
    real(dp) :: FDD, maxarg, g

    if (present(g_in)) then
       g = g_in
    else
       g = 1
    end if

    maxarg = 300. ! Upper cutoff for exponential function
    FDD = g/(1+exp(min((E-a)/b, maxarg)))
  end function FDD
    
    
  function reintegrand(p, k, w, kbT, mu, singulartol)
    ! The integrand in the integral to calculate the real part of the dielectric
    ! function. Assuming p, k, w are positive.
    ! 
    ! Parameters:
    ! ___________
    ! p: scalar
    !   The integration variable, which is also the momemtum of the electronic
    !   state.
    ! k: scalar
    !   The change of momentum for an incident photon with momentum k0 
    !   scattering to a state with momentum k1: k = |k1-k0|, in a.u.
    ! w: scalar
    !   The change of energy for an incident photon with energy w0 
    !   scattering to a state with energy w1: w = w0-w1, in a.u.
    ! kbT: scalar
    !   Thermal energy (kb - Boltzmann's constant, T is temperature) in a.u.
    ! mu: scalar
    !   Chemical potential in a.u.
    ! singulartol: scalar
    !   tolerance to handle singularities, where the arguments in the
    !   logarithms become 0.
    
    real(dp), intent(in) :: p, k, w, mu, kbT, singulartol
    real(dp) :: logpart, reintegrand, logtol, a, b, c, d, p1, p2

    logtol = 1e-3
    
    a = k**2 + 2*p*k + 2*w
    b = k**2 - 2*p*k + 2*w
    c = k**2 + 2*p*k - 2*w
    d = k**2 - 2*p*k - 2*w

    p1 = k**2 / (2*p*k + 2*w)
    p2 = k**2 / (2*p*k - 2*w)

    ! low "x" approximation for log(1+x)
    if (abs(p2) < logtol) then ! p1 < p2
       logpart = 2*p1 + 2*p2
    else
       ! tol avoids the singularities
       logpart = 0.5 * log((a**2 + singulartol)/(d**2 + singulartol)) &
            + 0.5*log((c**2 + singulartol)/(b**2 + singulartol))
       ! ! Equivalent expression
       ! logpart = log(abs(1+p1)) - log(abs(1-p1)) &
       !           + log(abs(1+p2)) - log(abs(1-p2))
       
    end if

    reintegrand = FDD(p**2/2, mu, kbT)*p*logpart    
    
  end function reintegrand

  function redielectric(k, w, kbT, mu, usertol, printmssg_in) result(er)
    ! Returns the real part of the RPA dielectric function.
    ! Integrates reintegrand while avoiding the singularities.
    !
    ! Parameters:
    ! ___________
    ! usertol: scalar
    !   A user set tolerance used to decide if an approximation of the
    !   integral is okay. This is to make sure our treatment of the singular
    !   points is fine.
    
    real(dp), intent(in) :: k, w, mu, kbT
    real(dp), intent(in), optional :: usertol
    logical, optional  :: printmssg_in ! Print warning message if .true.
    real(dp) :: xmax, x, dx, er, er2, tol, singular_tol, error
    integer  :: i, imax, loopcount
    logical :: printmssg

    if (present(usertol)) then
       tol = usertol
    else
       tol = 1e-4
    end if

    if (.not. present(printmssg_in)) then
       printmssg = .false.
    else
       printmssg = .true.
       end if
    
    singular_tol = 1e-6 ! singularity tolerance
    
    xmax = sqrt(2*(mu + 20*kbT)) ! 1/exp(20) \approx 1e-9
    imax = 2000
    dx = xmax/imax

    ! integrate
    er = 0.0
    do i=1,imax
       x = i*dx
       er = er + reintegrand(x, k, w, kbT, mu, singular_tol)*dx
    end do
    er2 = er ! treat er2 as a temporary place holder
    er = er2 + tol + 1 ! Need to get into do-while loop
    
    ! integrate at least twice to make sure our treatment of singularities is
    ! satisfactory.
    error = abs(er - er2)/er
    loopcount = 0
    do while ((er /= 0) .and. (error > tol)) ! Just a while loop
       ! Don't get stuck in a loop forever.
       if (loopcount >= 10) then
          ! We failed. Go home!
          if (printmssg) then
             print *, "RPAdielectric: Convergence failed due to singularities in ",&
                  "redielectric(): "
             print "(a, e9.4, a, e9.4, a, e9.2)", "k = ", k, "; w = ", w,&
                  "; confidence in result: ", error
          end if
          ! Return our best guess
          er = 1 + 2 / pi / k**3 * er2 
          return 
       end if
       loopcount = loopcount + 1
       
       er = er2 ! Hold previous er2 value
       er2 = 0.
       ! Try a smaller number for the singularities
       singular_tol = singular_tol/10
       
       do i=1,imax
          x = i*dx
          er2 = er2 + reintegrand(x, k, w, kbT, mu, singular_tol)*dx
       end do
       error = abs(er - er2)/er
    end do
    
    er = 1 + 2 / pi / k**3 * er2

  end function redielectric

  function imdielectric(k, w, kbT, mu) result(ei)
    ! Returns the imaginary part of the RPA dielectric.

    real(dp), intent(in) :: k, w, mu, kbT
    real(dp) :: ei, a, b

    a=abs(2*w-k**2)/2./k
    b=abs(2*w+k**2)/2./k
    ei=2 * kbT/k**3 * log(FDD(mu, b**2/2, kbT)/FDD(mu, a**2/2, kbT)) 
  end function imdielectric


  function elf(k, w, kbT, mu, printmssg)
    real(dp), intent(in) :: k, w, mu, kbT
    logical, intent(in), optional :: printmssg
    real(dp) :: elf, er, ei

    ! I'm surprised passing an optional arg unchecked worked (for now)...
    er = redielectric(k, w, kbT, mu, printmssg_in=printmssg)
    ei = imdielectric(k, w, kbT, mu)
    elf = ei / (ei**2 + er**2)
  end function elf
  
end module RPAdielectric
