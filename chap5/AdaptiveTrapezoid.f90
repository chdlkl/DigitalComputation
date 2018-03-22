Module Trapezoid_mod
  Implicit none
  Real(kind=8) :: a = 0.d0, b = 1.d0  !// 积分区间
  Real(kind=8) :: pi = acos(-1.d0), eps = 1.d-8, s0
  Real(kind=8) :: calPi = 0.d0  !// 数值积分解 
Contains
  Real(kind=8) function func ( x )
    Implicit none
    Real(kind=8) :: x
    
    func = 4.d0 / ( 1.d0 + x * x )  !// 此函数的积分值为pi
  End function func
  
  Recursive Real(kind=8) function trapezoid ( atmp, btmp, eps, s ) result( calPi )
    Implicit none
    Real(kind=8), intent(in) :: atmp, btmp, eps, s
    Real(kind=8) :: h, s1, s2
    
    h = ( btmp - atmp ) / 2.d0
    s1 = h * ( func(atmp) + func(atmp+h) ) / 2.d0
    s2 = h * ( func(atmp+h) + func(btmp) ) / 2.d0
    if ( abs( s - s1 - s2 ) < 3.d0*eps ) then
      calPi = s1 + s2
    else
      calPi = trapezoid( atmp, (atmp+btmp)/2.d0, eps/2.d0, s1 ) + trapezoid( (atmp+btmp)/2.d0, btmp, eps/2.d0, s2 )
    end if
  End function trapezoid
End module Trapezoid_mod
  
Program AdaptiveTrapezoid
  use Trapezoid_mod
  Implicit none
  !// 第一次计算trapezoid积分
  s0 = ( b-a ) * ( func(a) + func(b) ) / 2.d0 
  calPi = trapezoid( a, b, eps, s0 )
  Write( *, '(1x,a,2x,g0)') "the integral is", calPi
  Write( *, '(1x,a,e)') "the absolute error is", abs(pi-calPi)
End program AdaptiveTrapezoid  