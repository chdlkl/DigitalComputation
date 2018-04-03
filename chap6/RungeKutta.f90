!//---------------------------------------
!//本程序以四阶Runge-Kutta算法
!//u' = 4*t*sqrt(u), 0 =< t <= 2, u(0) = 1 精确解为：u(t) = (1+t^2)^2
!//---------------------------------------
Module Runge_Kutta
  Implicit none
  Real(kind=8), parameter :: t0 = 0.d0, t1 = 1.d0
  Real(kind=8), parameter :: h = 1d-2, num = nint((t1-t0)/h), u0=1.d0
  Real(kind=8) :: k1, k2, k3, k4
  Real(kind=8) :: t(0:Num),u(0:Num)
  Integer :: i
Contains
Subroutine solve()
  Implicit none
  u(0) = u0; t(0) = t0
  Do i = 1, num
    k1 = 2.d0 * t(i-1) * u(i-1)
    k2 = 2.d0 * ( t(i-1) + 0.5*h ) * ( u(i-1) + 0.5*h*k1 )
    k3 = 2.d0 * ( t(i-1) + 0.5*h ) * ( u(i-1) + 0.5*h*k2 )
    k4 = 2.d0 * ( t(i-1) + h ) * ( u(i-1) + h*k3 )
    t(i) = t0 + dble(i)*h
    u(i) = u(i-1) + h / 6.d0 * ( k1 + 2.d0*k2 + 2.d0*k3 + k4 )
  End do

  open( 100, file = 'shujv.txt', status='unknown' )
  Do i = 0, num
    write( 100,"(4f14.8)" ) t(i), exp( t(i)*t(i) ), u(i), abs( u(i)-exp(t(i)*t(i)) ) / exp( t(i)*t(i) ) * 100
  End do
  close(100)
End subroutine solve
End module Runge_kutta
  
Program main
  use Runge_Kutta  
  Implicit none
  call solve()
End program main