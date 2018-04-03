!// y' = t*y + t^3
!// y(t) = 3*exp(t^2/2) - t^2 - 2
Module Runge_Kutta
  Implicit none
  Real(kind=8), parameter :: t0 = 0.d0, t1 = 1.d0
  Real(kind=8), parameter :: h = 1d-2, num = nint((t1-t0)/h), u0 = 1.d0
  Real(kind=8) :: k1, k2, k3, k4
  Real(kind=8) :: t(0:num), u(0:num), u_tmp
  Integer :: i
Contains
Subroutine solve()
  Implicit none
  u(0) = u0; t(0) = t0
  Do i = 1, num
    k1 = derfunc( t(i-1), u(i-1) )
    k2 = derfunc( t(i-1)+h/2.d0, u(i-1)+h*k1/2.d0 )
    k3 = derfunc( t(i-1)+h/2.d0, u(i-1)+h*k2/2.d0 )
    k4 = derfunc( t(i-1)+h, u(i-1)+h*k3 )
    t(i) = t0 + dble(i)*h
    u(i) = u(i-1) + h / 6.d0 * ( k1 + 2.d0*k2 + 2.d0*k3 + k4 )
  End do

  open( 100, file = 'shujv.txt', status='unknown' )
  Do i = 0, num
    u_tmp = 3.d0 * exp(t(i)**2/2.d0) - t(i)**2 - 2.d0
    write( 100,"(4f14.8)" ) t(i), u_tmp, u(i), abs( u(i)-u_tmp ) / u_tmp * 100
  End do
  close(100)
End subroutine solve

Real(kind=8) function derfunc( t, u )
  Implicit none
  Real(kind=8), intent(in) :: t, u
  derfunc = t * u + t**3
End function derfunc
End module Runge_kutta
  
Program main
  use Runge_Kutta  
  Implicit none
  call solve()
End program main