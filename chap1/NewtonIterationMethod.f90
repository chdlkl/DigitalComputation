Program NewtonIterationMethod
  Implicit none
  Real(kind=8), external :: func, DerivativeFunc
  Real(kind=8), parameter :: eps = 1.d-8 
  Real(kind=8) :: x0, xtmp
  Integer, parameter :: nloop = 1000 !// 迭代次数上限
  Integer :: i 
  
  i = 0
  x0 = -0.7  !// 初始估计值
  Do 
    xtmp = x0 - func(x0) / DerivativeFunc(x0)        !// X(i+1) = Xi - f(Xi) / f'(Xi)  i = 1,2,...N
    i = i + 1
    If ( abs( xtmp-x0 ) < eps .and. i > nloop ) exit
    x0 = xtmp
  End do 
  Write ( *,'(a,g0)' ) ' x^3 + x - 1的零点为:', xtmp
  
End program NewtonIterationMethod
    
Real(kind=8) function func( x0 )  !// 原函数：f(x) = x^3 + x - 1 
  Implicit none 
  Real(kind=8), intent(in) :: x0 
  func = x0**3 + x0 - 1.d0
End function func 
  
Real(kind=8) function DerivativeFunc( x0 ) result( res )  !// 原函数导数：f'(x) = 3x^2 + 1
  Implicit none 
  Real(kind=8), intent(in) :: x0 
  res = 3.d0 * x0**2 + 1.d0
End function DerivativeFunc