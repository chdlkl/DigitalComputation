Program FixedPointIteration  !// 使用不动点迭代方法求解方程 x^3 + x - 1 = 0 的零点
  Implicit none 
  Real(kind=8), parameter :: eps = 1.d-8
  Real(kind=8), external :: func
  Real(kind=8) :: xtemp, x0 = 0.5d0  !// x0为迭代初始值
  
  Do 
    xtemp = func ( x0 )
    If ( abs( xtemp - x0 ) < eps ) exit 
    x0 = xtemp
  End do 
  Write ( *,'(a,f13.9)' ) ' the result of fixed point iteration is', xtemp

End program FixedPointIteration
  
Real(kind=8) function func ( x )
  Implicit none 
  Real(kind=8), intent(in) :: x 
  func = ( 1.d0 + 2.d0 * x * x * x ) / ( 1.d0 + 3.d0 * x * x )
End function func 
!// 方程 x^3 + x - 1 = 0 主要的不动点格式有三种
!// 第一种：x = 1 - x^3 这一种格式不收敛
!// 第二种：x = ( 1 - x^3 )^(1/3)，这种格式收敛
!// 第三种：就是此代码中的收敛格式。x = ( 1 + 2*x^3 ) / ( 1 + 3*x^2 )，这种格式不仅收敛，而且收敛速度快，精度高