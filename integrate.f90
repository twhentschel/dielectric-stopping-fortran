module integrate

  use utils, only : dp
  
  implicit none
  
contains

  pure function trapezoidal(y, x) result(r)
    ! Code taken from Fortran Wiki website:
    ! http://fortranwiki.org/fortran/show/integration
    ! written by contributor jabirali on April 22, 2016.
    !
    ! Calculates the integral of an array y with respect to x using the
    ! trapezoidal approximation. Note that the mesh spacing of x does not have
    ! to be uniform.
    real(dp), intent(in)  :: x(:)         !! Variable x
    real(dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(dp)              :: r            !! Integral ∫y(x)·dx
    integer               :: n

    n = size(x)
    ! Integrate using the trapezoidal rule
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    ! associate(n => size(x))
    !   r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    ! end associate
  end function trapezoidal
  
end module integrate


! program main
!   use integrate
!   implicit none

!   real(dp) :: f(100), x(size(f)), a, b
!   integer  :: i, n

!   n = size(f)
!   ! integrate x^2 from 0 to 1
!   a = 0.
!   b = 1.
!   do i = 1,n
!      x(i) = (i-1) * (b-a) / (n-1) + a
!   end do
!   f(1:n) = x(1:n)**2

!   print "(f7.4)", trapezoidal(f, x)
  
! end program main
