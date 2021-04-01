module RPAdielectric
  implicit none
  private
  
  integer, parameter :: dp = selected_real_kind(15, 307)
  
contains
  
  function FDD(E, a, b)
    ! Fermi-Dirac Distriution function.
    real(dp), intent(in) :: E, a, b
    real(dp) :: FDD

    FDD = 1/(1+exp((E-a)/b))
  end function FDD
    
    
  function reintegrand(p, k, w, mu, kbT, tol)
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
    ! tol: scalar
    !   Some user defined tolerance
    
    real(dp), intent(in) :: p, k, w, mu, kbT, tol
    real(dp) :: logpart, reintegrand
    real(dp) :: a, b, c, d, p1, p2
    ! interface
    !    function FDD(E, mu, kbT)
    !      import
    !      real(dp), intent(in) :: E, mu, kbT
    !      real(dp) :: FDD
    !    end function FDD
    ! end interface
    
    a = k**2 + 2*p*k + 2*w
    b = k**2 - 2*p*k + 2*w
    c = k**2 + 2*p*k - 2*w
    d = k**2 - 2*p*k - 2*w

    p1 = k**2 / (2*p*k + 2*w)
    p2 = k**2 / (2*p*k - 2*w)

    ! low "x" approximation for log(1+x)
    if (p2 < tol) then ! p1 < p2
       logpart = -2*p1 - 2*p2
    else
       ! logpart = log(abs(1-p1)) - log(abs(1+p1)) &
       !           + log(abs(1-p2)) - log(abs(1+p2))
       logpart = log(abs(a*c/(b*d))
    end if


    
    reintegrand = FDD(p**2/2, mu, kbT)*p*logpart    
    
  end function integrand
