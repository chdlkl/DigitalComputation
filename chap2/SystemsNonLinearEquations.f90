!// 方程组：-u^3 + v = 0
!//         u^2 + v^2 - 1 = 0
!// 即：f1(u,v) = -u^3 + v  f2(u,v) = u^2 + v^2 - 1
!// 多元牛顿法的求解适用于小型方程组，对于大型不适用，因为需要提前指导导数
!// 多变量牛顿法
!// x0 = 初始向量
!// DF(Xk) * s = -F(Xk)
!// Xk+1 = Xk + s, k = 0, 1, 2, ...
!// 对于多元牛顿法来说，当方程组的阶数较小时，求解DF与F相对容易
!// 但当方程组的阶数过大时，求解DF与F相当繁琐。
Program SystemsNonLinearEquations
  Implicit none
  Integer :: i
  Integer, parameter :: m = 2, maxLoop = 50
  Real(kind=8), parameter :: eps = 1.d-12
  Real(kind=8) :: DF(m,m) = 0.d0, F(m,1) = 0.d0, s(m,1) = 0.d0
  Real(kind=8) :: x(m,1) = [ 1.d0, 2.d0 ], x0(m,1)  !// 赋初值
  
  x0 = x
  Do i = 1, maxLoop
    call GetDF ( DF, m, x0 )
    call GetF ( F, m, x0 )
    call GaussianElimination (  DF, F, s, m )
    x = x0 + s
    If ( maxval(abs(x-x0)) < eps ) exit
    x0 = x
  End do
  
  Write ( *,'(1x,a,g0)' ) '迭代次数为：', i
  Write ( *,'(1x,a)' ) '非线性方程组的解为：'
  Do i = 1, m
    Write ( *,'(g0)' ) x(i,1)
  End do
End program SystemsNonLinearEquations
  
Subroutine GetDF ( DF, m, x )
  Implicit none
  Integer, intent(in) :: m
  Real(kind=8), intent(in) :: x(m,1)
  Real(kind=8), intent(inout) :: DF(m,m)
  Real(kind=8) :: coff(m,m), coffuv(m,m)
  
  coff = reshape( [ -3.d0, 2.d0, 1.d0, 2.d0 ], [m,m] )
  coffuv(1,1) = x(1,1)**2; coffuv(2,1) = x(1,1)
  coffuv(1,2) = 1.d0; coffuv(2,2) = x(2,1)
  DF = coff * coffuv

End subroutine GetDF
  
Subroutine GetF ( F, m, x )
  Implicit none
  Integer, intent(in) :: m
  Real(kind=8), intent(inout) :: F(m,1)
  Real(kind=8), intent(in) :: x(m,1)
  Real(kind=8) :: f1, f2
  
  !// f1 = -u^3 + v
  !// f2 = u^2 + v^2 - 1
  F(1,1) = -1.d0 * ( -x(1,1)**3 + x(2,1) )
  F(2,1) = -1.d0 * ( x(1,1)**2 + x(2,1)**2 - 1.d0 )
End subroutine GetF
  
Subroutine GaussianElimination ( a, b, x, m )  !// 高斯消去
  Implicit none 
  Integer :: i, j, k 
  Real(kind=8), parameter :: eps = 1.d-4  !// 当主元小于这个数时，程序退出，(可用部分主元法进行改善!)
  Real(kind=8) :: mult
  Integer, intent(in) :: m
  Real(kind=8), intent(inout) :: a(m,m), b(m,1), x(m,1)
  
  Do j = 1, m - 1
    If ( abs(a(j,j)) < eps ) Then
      Write ( *,'(1x,a)' ) ' The pivot is zero!'
      stop 
    End if 
    Do i = j + 1, m
      mult = a(i,j) / a(j,j)
      Do k = j, m 
        a(i,k) = a(i,k) - mult * a(j,k)
      End do 
      b(i,1) = b(i,1) - mult * b(j,1)
    End do 
  End do  
  !// 回代
  Do i = m, 1, -1
    Do j = i + 1, m 
      b(i,1) = b(i,1) - a(i,j) * x(j,1)
    End do 
    x(i,1) = b(i,1) / a(i,i)
  End do 
  
End subroutine GaussianElimination 
