module MerminDielectric
  use utils, only : dp, pi
  use RPAdielectric, only : redielectric, imdielectric
  implicit none
  private
  public mrpac, mELF

contains

  function mELF(k, w, kbT, mu, vr, vi) result(elf)
    real(dp), intent(in) :: k, w, kbT, mu, vr, vi
    real(dp) :: elf, er, ei
    call mrpac(k, w, kbT, mu, vr, vi, er, ei)
    elf = ei / (er**2 + ei**2)
  end function mELF
  
  subroutine mrpac(k,w,t,m,vr,vi,er,ei)
    !returns real and imaginary part of Mermin RPA dielectric er and ei following T Hentschel code
    !for k,w,t, m, and v in atomic units through direct numerical integration
    !input k is momentum, w is energy, t is temperature, m is chemical potential
    !vr and vi are real and imaginary parts of the dynamic collision frequency at w 
    real(dp), intent(in) ::  k,w,t,m
    real(dp), intent(out) :: er, ei
    real(dp) :: p, dpp,x,pmax,f,a,b,c,d,vr,vi,phi,del, xmx
    complex ::  v,em,e0,ev,ci
    integer  :: i,imax

    del = 1e-8
    imax = 2000
    xmx = 300
    pmax=2.*max(m+3.*t,5.*t) !energy max
    pmax=sqrt(2.*pmax)       !momentum max
    dpp=pmax/real(imax)
    !integrate real & imaginary part of dielectric
    er=0.0
    ei=0.0
    do i=1,imax
       p=real(i)*dpp
       x=(p*p/2.-m)/t
       f=1.0/(1.0+exp(min(xmx,x)))
       x=2.*vr
       a=k*k+2.*p*k+2.*(w-vi)
       b=k*k-2.*p*k+2.*(w-vi)
       c=k*k+2.*p*k-2.*(w-vi)
       d=k*k-2.*p*k-2.*(w-vi)
       phi= atan2(x,a) - atan2(x,b) + atan2(-x,c) - atan2(-x,d)
       ei=ei+f*p*phi*dpp
       x=x*x
       a=a*a+x
       b=b*b+x
       c=c*c+x
       d=d*d+x
       a=sqrt(abs(a/b))
       b=sqrt(abs(c/d))
       er=er+f*p*(log(a)+log(b))*dpp
    end do
    er=1.+2./pi/k**3*er
    ei=2./pi/k**3*ei
    
    
    
    !Mermin (complex!!)
    ev=cmplx(er,ei, kind=dp) !real and imaginary parts of e(k, w + iv)
    er = redielectric(k, 0._dp, t, m)
    ei = imdielectric(k, 0._dp, t, m)
    e0=cmplx(er,ei, kind=dp)   !real and imaginary of e(k,0)
    v=cmplx(vr,vi, kind=dp)
    ci=cmplx(0.,1., kind=dp)
    em= 1. + (w+del+ci*v)*(ev-1.)/(w +del + ci*v*(ev-1.)/(e0-1.)) !Tommy
    er=real(em)
    ei=aimag(em)
    return
    
  end subroutine mrpac
  
end module MerminDielectric
