module ArrayOps

  use utils, only : dp
  
  implicit none

contains

  function linspace(start, end, N) result(samples)
    real(dp), intent(in)  :: start, end
    integer, intent(in)   :: N
    real(dp)              :: samples(N), step
    integer               :: i

    step = (end - start) / real(N - 1, dp)
    do i = 1,N
       samples(i) = (i-1) * step + start
    end do
  end function linspace

  function logspace(start, end, N) result(samples)
    real(dp), intent(in) :: start, end
    integer , intent(in) :: N
    real(dp)             :: samples(N)
    
    samples = exp(linspace(log(start), log(end), N))
    
  end function logspace
  
  function bunchspace(center, width, N) result(samples)
    ! Like logspace, but the bunch is centered around center, not 0.
    ! It is symmetric and extends from [center - width, center + width], and
    ! puts this in arr. Total number of points will be 2N-1.
    ! center must be greater than width.
    real(dp), intent(in)    :: center, width
    integer,  intent(in)    :: N
    real(dp)                :: logsamples(N), samples(2*N-1), begin

    begin = center - width

    ! Treat first part of samples as a temporary holder at first
    logsamples = exp(linspace(log(begin), log(center), N))
    samples(N:1:-1) = (center + begin) - logsamples
    samples(N+1:)= width + logsamples(2:)

  end function bunchspace

  function binarysearch(arr, val) result(index)
    ! Search through a real array arr to look for the value val using the binary
    ! search algorithm. If val does not exist in arr, return the index of the
    ! element in arr that occurs immediately AFTER val.

    real(dp), intent(in) :: arr(:), val
    integer              :: index, L, R, mid

    L = 1
    R = size(arr)
    do while ((L <= R))
       mid = floor((L+R)/ real(2,dp))
       if (arr(mid) < val) then
          L = mid + 1
       else if (arr(mid) > val) then
          R = mid - 1
       else
          index = mid
          return
       end if
    end do
    index = L
    return
  end function binarysearch

  ! function sortedinsert(arr, val) result(newarr)
  !   ! Insert a value (val) into the sorted array (arr).
  !   ! Array is assumed to sorted from least to greatest.
  !   real(dp), intent(in) :: arr(:), val
  !   real(dp)             :: newarr(size(arr) + 1)
  !   integer              :: N, i

    
    
  ! end function sortedinsert
  
end module ArrayOps


! program main
!   use ArrayOps

!   implicit none

!   real(dp) :: A(11)
!   integer  :: i
  
!   A = logspace(0.01_dp, 10._dp, size(A))
!   print '(f7.4)', A
!   !i = binarysearch(A, -1._dp)
!   !print *, i


  
! end program main
