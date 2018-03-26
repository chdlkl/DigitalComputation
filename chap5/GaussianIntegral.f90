!// 注: 参考来源http://fcode.cn/algorithm-73-1.html  
Module Gauss_Legendre !//高斯—勒让德积分高斯点及权重的求解模块
  Implicit none
  Integer, parameter :: n  = 11                          !// 设置求解高斯点的个数
  Integer, parameter :: DP = selected_real_kind( p=13 )  !// 设置kind的数值
  Real(kind=DP), parameter :: eps = 1.0e-15_DP           !// 精度设置
Contains
  Real(Kind=DP) function N_Legendre(x) !// 生成n阶勒让德多项式
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a(n), x
    a(1) = x !// 1阶勒让德多项式
    a(2) = 1.5_DP*x*x - 0.5_DP !// 2阶勒让德多项式
    Do i = 3, n
      a(i) = ( dble(i+i-1)*x*a(i-1) - dble(i-1)*a(i-2) ) / dble(i) !// 利用递推关系产生n阶勒让德多项式
    End do
    N_Legendre=a(n) !//生成的n阶勒让德多项式
  End function N_Legendre

  Real(Kind=DP) Function N1_Legendre(x)  !// 生成n-1阶勒让德多项式 
    Implicit none
    Integer :: i
    Real (Kind=DP) :: a(n), x
    a(1) = x
    a(2) = 1.5_DP*x**2 - 0.5_DP
    Do i = 3, n - 1
      a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
    End Do
    N1_Legendre = a(n-1)     
  End function N1_Legendre
  
  Real(Kind=DP) function DN_Legendre(x)  !// 生成n阶勒让德多项式的导数表达式
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a(n), x
    a(1) = x  !// 1阶勒让德多项式
    a(2) = 1.5_DP*x*x - 0.5_DP !// 2阶勒让德多项式
    Do i = 3, n
      a(i) = ( dble(i+i-1)*x*a(i-1) - dble(i-1)*a(i-2) ) / dble(i) !// 利用递推关系产生n阶勒让德多项式
    End Do
    DN_Legendre = ( a(n-1) - x*a(n) )*dble(n) / (1.0_DP - x*x ) 
  End function DN_Legendre

  Real(Kind=DP) function NewtonIteration(a, b) !// 牛顿法求解函数的解
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a, b, x, xtmp
    Integer, parameter :: nloop = 2000
    !// a,b是传递进来的划分好的有一个解存在的区间
    x = ( a + b ) / 2.d0  !// 初始估计值
    i = 0
    Do 
      xtmp = x - N_Legendre(x) / DN_Legendre(x)   !// X(i+1) = Xi - f(Xi) / f'(Xi)  i = 1,2,...N
      i = i + 1
      If ( abs( xtmp-x ) < eps .and. i > nloop ) exit
      x = xtmp
    End do 
    NewtonIteration = x
  End function NewtonIteration

  Subroutine root_coeff ( f_root, f_coeff )  !// 计算N阶勒让德多项式的根与去做权重系数
    Implicit none
    Real(Kind=DP) :: m, nstep, f_root(n), f_coeff(n) !// 定义数组,大小n由module开始声明。
    Integer :: i, j
    Real(kind=DP), parameter :: h = 1.d-6
    j = 0   !// 赋值控制循环变量的初值           
    m = -1.d0 + h   !// 设置计算域[-1，1] 的下限，即代替-1 
    nstep = nint(2.d0/h)
    Do i = 1, nstep   !// 这个循环次数应该是由步长0.000001决 定,计算方法：2000000=2/0.000001     
      If ( N_Legendre(m)*N_Legendre(m+h) < 0 ) then   !// 从下限处开始往上逐步累加
        j = j + 1    !// 记录这是第几个解
        f_root(j) = NewtonIteration( m, m+h )!// 调用牛顿法求解程序在分好的一小段上求解，将解存储在fn（j）
        f_coeff(j) = 2.0_DP / ( dble(n) * N1_Legendre(f_root(j)) * DN_Legendre(f_root(j)) ) !// 利用公式计算高斯点的权重
        write (*,'(1x,a,g0)') '高斯点序号: ', j
        write (*,'(1x,a,g0,2x,a,g0)') '高斯点: ', f_root(j), '高斯点权重', f_coeff(j)
        write(*,'(1x,a)') '------------------------------------------------------'
      End if
      m = m + h !// 执行完一次判断m向前推进一步
    End Do
  End subroutine root_coeff
  
  Real(Kind=DP) function func(x) !// 被积函数
    Implicit none
    Real(Kind=DP) :: x
    func = exp( -x*x/2.0_DP ) !// 每次计算在这里修改被积函数即可
  End function func

End module Gauss_Legendre

Program GaussianIntegral
  use Gauss_Legendre
  Implicit none
  Real (Kind=DP) :: f_root(n), f_coeff(n), x, a, b, answer
  Integer :: i
  Call root_coeff ( f_root, f_coeff ) !// 调用求高斯零点和权重的子函数
  a = -1.0_DP !// 积分上限
  b = 1.0_DP !// 积分下限   
  answer = 0.d0 !// 求积分结果赋初始值
  !// 一般区间[a,b]上的积分公式
  !// Integral[f(x),a,b] = Integral[f( ((b-a)*t+b+a)/2 )] * ( b-a )/2. t为N阶勒让德多项式的根
  Do i = 1, n
    answer = answer + f_coeff(i) * func( (a+b) / 2.0_DP + (b-a) / 2.0_DP * f_root(i) ) !// 高斯勒让德求积分公式     
  End Do
  answer = answer * (b-a) / 2.0_DP
  !// 精确解为1.71124878378430
  !// 数值解为1.711248783784299  
  Write(*,'(1x,a,g0)') '高斯-勒让德求积分结果: ', answer

End program GaussianIntegral
