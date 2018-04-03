!// 此代码为RungeKuttaFehlberg的半成品
!// y' = t*y + t^3
!// y(t) = 3*exp(t^2/2) - t^2 - 2
Module RungeKuttaFehlberg
  Implicit none
  Real(kind=8), parameter :: t0 = 0.d0, t1 = 1.d0
  Real(kind=8), parameter ::  u0 = 1.d0
  Real(kind=8) :: k1, k2, k3, k4, k5, h = 1d-1  !// 给定初始步长
  Real(kind=8) :: t, u, u_tmp
  Integer :: i
Contains
Subroutine solve()
  Implicit none
  u = u0; t = t0
  open( 100, file = 'shujv.txt', status='unknown' )
  write( 100,"(4f14.8)" ) t, u, u0, abs( u-u0 ) / u0 * 100
  Do 
    k1 = derfunc( t, u )
    k2 = derfunc( t + h/4.d0, u + h*k1/4.d0 )
    k3 = derfunc( t + 3.d0*h/8.d0, u + 3.d0*h*k1/32.d0 + 9.d0*h*k2/32.d0 )
    k4 = derfunc( t + 12.d0*h/13.d0, u + 1932.d0*h*k1/2197.d0 - 7200.d0*h*k2/2197.d0 + 7296.d0*h*k3/2197.d0 )
    k5 = derfunc( t + h, u + 439.d0*h*k1/216.d0 - 8.d0*h*k2 + 3680.d0*h*k3/513.d0 - 845.d0*h*k4/4104.d0 )
    u = u + h * ( 25.d0*k1/216.d0 + 1408.d0*k3/2565.d0 + 2197.d0*k4/4104.d0 - k5/5.d0 )
    t = t + h
    If ( t > t1 ) exit
    u_tmp = 3.d0 * exp(t**2/2.d0) - t**2 - 2.d0
    write( 100,"(4f14.8)" ) t, u, u_tmp, abs( u-u_tmp ) / u_tmp * 100
  End do
  close(100)
End subroutine solve

Real(kind=8) function derfunc( t, u )
  Implicit none
  Real(kind=8), intent(in) :: t, u
  derfunc = t * u + t**3
End function derfunc
End module RungeKuttaFehlberg
  
Program main
  use RungeKuttaFehlberg  
  Implicit none
  call solve()
End program main