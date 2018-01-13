Program SecantMethod
  Implicit none 
  Real(kind=8) :: x0 = 0.d0, x1 = 1.d0
  Real(kind=8) :: xtmp
  Real(kind=8), parameter :: eps = 1.d-8
  Real(kind=8), external :: func
  Integer, parameter :: nloop = 1000 !// 设置最大循环次数
  Integer :: i

  i = 0
  Do
    If ( i > nloop ) exit
    xtmp = x1 - func(x1) * ( x1 - x0 ) / ( func(x1) - func(x0) ) !// Xi+1 = Xi - f(Xi) * ( Xi - Xi-1 ) / ( f(Xi) - f(Xi-1) ), i = 1, 2, 3, 4......
    If ( abs(xtmp - x1) <= eps ) exit 
    x0 = x1
    x1 = xtmp
    i = i + 1
  End do 
  Write ( *,'(a,g0)' ) ' The root of x^3 + x - 1 is ', xtmp

End Program SecantMethod

Real(kind=8) function func (x) !// 求方程 x^3 + x - 1 = 0 的零点
  Implicit none 
  Real(kind=8), intent(in) :: x 
  func = x**3 + x - 1.d0 
End function
