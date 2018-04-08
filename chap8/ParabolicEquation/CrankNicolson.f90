!// Ut = DUxx   a<=x<=b, t>=0  本例中D=4
!// U(x,0) = f(x)  a<=x<=b
!// U(a,t) = L(t)  t>=0
!// U(b,t) = R(t)  t>=0
Module BckDeff
  Implicit none
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0, D = 4.d0
  Real(kind=8), parameter :: y0 = 0.d0, yT = 1.d0
  Real(kind=8), parameter :: pi = acos(-1.d0)
  Integer, parameter :: numH = 10, numT = 100 !// 空间步数与时间步数
Contains
  Subroutine solve()
    !Include 'link_fnl_shared.h'  !// 安装imsl函数库的可以调用求解逆矩阵, 去掉13，14，56，65，66行的注释，并注释67行即可
    !use LINRG_INT
    Implicit none
    Integer :: i, j
    Real(kind=8) :: h, k, x, t, sigma  !// D为扩散系数，h为空间步长，k为时间步长
    Real(kind=8) :: a(numH-1,numH-1), b(numH-1,numH-1), u(0:numH,0:numT), Lside(0:numT), Rside(0:numT)
    Real(kind=8) :: tmp_u(numH-1,1), slu(numH-1,1), tmp(numH-1,1)
    
    h = ( xR-xL ) / numH
    k = ( yT-y0 ) / numT
    sigma = D*k / ( h*h )
    
    !// 构造系数矩阵a
    a = 0.d0
    Do i = 1, numH-1
      a(i,i) = 2.d0 + 2.d0*sigma
    End do
    Do i = 1, numH-1 - 1
      a(i,i+1) = -sigma
      a(i+1,i) = -sigma
    End do
    
    !// 构造系数矩阵b
    b = 0.d0
    Do i = 1, numH-1
      b(i,i) = 2.d0 - 2.d0*sigma
    End do
    Do i = 1, numH-1 - 1
      b(i,i+1) = sigma
      b(i+1,i) = sigma
    End do
    
    !// 边界条件
    Do i = 0, numT
      t = y0 + dble(i) * k
      Lside(i) = exp(t)
      Rside(i) = exp(t-0.5d0)
    End do
    
    !// 初值f(x,0) = exp(-x/2)
    u(1:numH-1,0) = exp( -(xL + [1:numH-1]*h) / 2.d0 )
    
    !// 计算
    !call LINRG( a, Inva )
    Do j = 0, numT-1
      tmp(1,1) = Lside(j) + Lside(j+1)
      tmp(numH-1,1) = Rside(j) + Rside(j+1)
      tmp(2:numH-2,1) = 0.d0
      tmp_u(1:numH-1,1) = u(1:numH-1,j)
      tmp_u = matmul( b,tmp_u ) + sigma * tmp
      slu = 0.d0
      call gauss( a, slu, tmp_u, numH-1 )  !// 利用gauss消去法进行求解
      !tmp_u = matmul( Inva, tmp_u )
      !u(1:numH-1,j+1) = tmp_u(1:numH-1,1)
      u(1:numH-1,j+1) = slu(1:numH-1,1)
    End do
    u(0,:) = Lside
    u(numH,:) = Rside
    
    !// output
    open( 100, file = 'output.dat' )
    Do i = 0, numH
      x = xL + dble(i) * h
      Do j = 0, numT
        t = y0 + dble(j) * k
        write( 100,'(4f15.8)' ) x, t, u(i,j), exp( t - x/2.d0 )   !// 解析解u(x,t) = exp( t - x/2.d0 ) 
      End do
    End do
    close( 100 )
    
  End subroutine solve
  
  Subroutine gauss( a, y, b, n )
    Implicit none
    Integer, intent(in) :: n
    Real(kind=8), intent(in) :: a(n,n)
    Real(kind=8), intent(out) :: y(n,1), b(n,1)
    Integer :: i, j, k
    Real(kind=8) :: mult, aa(n,n)
    
    aa = a
    Do j = 1, n - 1
      Do i = j + 1, n
        mult = aa(i,j) / aa(j,j)
        Do k = j, n 
          aa(i,k) = aa(i,k) - mult * aa(j,k)
        End do 
        b(i,1) = b(i,1) - mult * b(j,1)
      End do 
    End do 
    
    !// 回代
    y = 0.d0
    Do i = n, 1, -1
      Do j = i + 1, n 
        b(i,1) = b(i,1) - aa(i,j) * y(j,1)
      End do 
      y(i,1) = b(i,1) / aa(i,i)
    End do  
  End subroutine gauss
  
End module BckDeff

Program main
  use BckDeff
  Implicit none
  call solve
End program main  