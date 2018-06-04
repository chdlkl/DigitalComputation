Module mod
  Implicit none
  Real(kind=8), parameter :: a = 1.d0, b = 1.02d0  !// 积分区间
  Real(kind=8), parameter :: eps = 1.d-15  !// 误差控制
  Real(kind=8) :: s  !// 最终积分结果
  Integer :: n  !// 最终分的份数
Contains
Subroutine solve( s, a, b, tol, n )
!//  Purpose   :  自动变步长Simpson积分方法函数
!    
!//  Input  parameters  :
!//       1.  func 外部函数
!//       2. 
!//       3.  a,b积分区间
!//       4.  tol  积分误差容限
!
!//  Output parameters  :
!//       1.   s  积分结果
!//       2.   n  实际区间划分个数
  Real(kind=8) :: s, s1, a, b, del, tol
  Integer :: n, i, m

  !//初始划分40个子区间
  n = 40  !//n可以根据被积函数的复杂程序适当地减小或增大

  !//最大允许划分100次
  m = 100  !//此处的m也可以根据被积函数的复杂程序适当地减小或增大
  Do i = 1, m
    call simp( s, a, b, n ) 
    n = n * 2
    call simp( s1, a, b, n )
    del = abs(s-s1)

    !//满足精度后就停止循环
    if ( del < tol ) exit
  End do

  s = s1

end subroutine solve

Subroutine simp( s, a, b, n )
  Real(kind=8) :: s, a, b, h, f1, f2, f3, f4, t1, t2
  Integer :: n, k

  s = 0d0
  h = (b-a) / n / 2d0

  call func1( f1, a )
  call func1( f2, b )

  s = f1 + f2

  !//k=0 情况
  call func1( f1, a+h )
  s = s + 4d0*f1

  Do k = 1, n-1
    t1 = a + (2d0*k+1.d0) * h
    t2 = a + 2d0 * k * h

    call func1( f3, t1 )
    call func1( f4, t2 )

    s = s + f3*4d0 + f4*2d0  
  End do

  s = s*h / 3d0

End subroutine simp


Subroutine func1( f, x )
  Implicit none
  Real(kind=8) :: f, x
  f = -10.d0 * x**(-11.d0) !//被积函数(原函数的导数)
End subroutine func1

End module mod

Program main  !// 此程序一般只需更改积分区间[a,b],误差限度eps,以及被积函数fun1
  use mod
  call solve( s, a, b, eps, n )
  write( *,'(1x,"区间划分等份为:",I5,/,1x,"积分结果为:",g0)' ) n, s
  write( *,'(1x,"精确积分为:",g0)' ) 1.02d0**(-10.d0) - 1.d0
End program main
