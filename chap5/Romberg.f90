!// author: luk
!// time: 2018/3/20
!// purpose: Romberg integral for ln(x) between 1 and 2
Module mod
  Implicit none
Contains
  Subroutine cal_integral( )
    Implicit none
    Integer :: i, j, k
    Integer, parameter :: n = 4  !// 步长减半的次数
    Real(kind=8) :: r(n,n), h
    Real(kind=8) :: a = 1.d0, b = 2.d0  !// 积分上下限
    Real(kind=8) :: total
    
    r = 0.d0
    r(1,1) = ( b - a ) * ( func(a) + func(b) ) / 2.d0
    Do j = 2, n
      total = 0.d0
      h = ( b - a ) / 2**( j - 1 )
      Do i = 1, 2**( j-2 )
        total = total + func( a + (2*i-1)*h )
      End do
      r(j,1) = r(j-1,1) / 2.d0 + h * total
      Do k = 2, j
        r(j,k) = ( 4.d0**(k-1) * r(j,k-1) - r(j-1,k-1) ) / ( 4.d0**(k-1) - 1.d0 )
      End do
    End do  
    
    !// output r
    Do i = 1, n
      Write(*,'(*(f20.14))') r(i,:)
    End do
    Write(*,'(1x,a)') "精确积分为:"
    Write(*,'(*(f20.14))') r(n,n)
    
  End subroutine cal_integral
  
  Real(kind=8) function func( x )
    Implicit none
    Real(kind=8) :: x
    
    func = log(x)
  End function func
End module mod
  
Program Romberg
  use mod
  Implicit none
  call cal_integral
End program Romberg