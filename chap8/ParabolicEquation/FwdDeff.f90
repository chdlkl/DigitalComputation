!// Ut = DUxx   a<=x<=b, t>=0
!// U(x,0) = f(x)  a<=x<=b
!// U(a,t) = L(t)  t>=0
!// U(b,t) = R(t)  t>=0
Module FwdDeff
  Implicit none
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0
  Real(kind=8), parameter :: y0 = 0.d0, yT = 1.d0
  Real(kind=8), parameter :: pi = acos(-1.d0)
  Integer, parameter :: numH = 10, numT = 250 !// 空间步数为10，因为两个端点的值已知，所以计算的点只有9个
Contains
  Subroutine solve()
    Implicit none
    Integer :: i, j
    Real(kind=8) :: D, h, k, x, t, sigma  !// D为扩散系数，h为空间步长，k为时间步长
    Real(kind=8) :: a(numH-1,numH-1), u(0:numH,0:numT), Lside(0:numT), Rside(0:numT)
    Real(kind=8) :: tmp_u(numH-1,1), b(numH-1,1)
    
    D = 1.d0
    h = ( xR-xL ) / numH
    k = ( yT-y0 ) / numT
    sigma = D*k / ( h*h )
    
    !// 构造系数矩阵a
    a = 0.d0
    Do i = 1, numH-1
      a(i,i) = 1.d0 - 2.d0*sigma
    End do
    Do i = 1, numH-1 - 1
      a(i,i+1) = sigma
      a(i+1,i) = sigma
    End do
    
    !// 边界条件
    Do i = 0, numT
      t = y0 + dble(i) * k
      Lside(i) = 0.d0 * t
      Rside(i) = 0.d0 * t
    End do
    
    !// 初值条件
    !// f(x,0) = sin(2*pi*x)**2
    u(1:numH-1,0) = sin( 2.d0*pi* ( xL + [1:numH-1]*h ) )**2
    
    !// 计算
    Do j = 0, numT-1
      b(1,1) = Lside(j)
      b(numH-1,1) = Rside(j)
      b(2:numH-2,1) = 0.d0
      tmp_u(1:numH-1,1) = u(1:numH-1,j)
      tmp_u = matmul( a,tmp_u ) + sigma * b
      u(1:numH-1,j+1) = tmp_u(1:numH-1,1)
    End do
    u(0,:) = Lside
    u(numH,:) = Rside
    
    !// output
    open( 100, file = 'output.dat' )
    Do i = 0, numT
      t = y0 + dble(i) * k
      Do j = 0, numH
        x = xL + dble(j) * h
        write( 100,'(3f15.8)' ) x, t, u(j,i)
      End do
    End do
    close( 100 )
    
  End subroutine solve
End module FwdDeff

Program main
  use FwdDeff
  Implicit none
  call solve
End program main  