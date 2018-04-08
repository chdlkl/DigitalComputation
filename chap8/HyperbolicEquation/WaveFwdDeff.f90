!// Utt = c^2*Uxx   a<=x<=b, t>=0  c = 2
!// U(x,0) = f(x) = sin(pi*x)  a<=x<=b
!// Ut(x,0) = g(x) = 0  a<=x<=b
!// U(a,t) = L(t)  t>=0
!// U(b,t) = R(t)  t>=0
Module FwdDeff
  Implicit none
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0, c = 2.d0
  Real(kind=8), parameter :: y0 = 0.d0, yT = 1.d0
  Real(kind=8), parameter :: pi = acos(-1.d0)
  Integer, parameter :: numH = 100, numT = 2000 !// 空间步数和时间步数
Contains
  Subroutine solve()
    Implicit none
    Integer :: i, j
    Real(kind=8) :: h, k, x, t, sigma  !// D为扩散系数，h为空间步长，k为时间步长.  h与k要满足: c*k <= h称为波方程的CFL条件
    Real(kind=8) :: u(0:numH,0:numT), Lside(0:numT), Rside(0:numT)
    
    h = ( xR-xL ) / numH
    k = ( yT-y0 ) / numT
    sigma = c*k / h
    If ( c*k >= h ) then
      write( *,'(1x,a)' ) ' 不满足CFL条件, 程序中止... '
      stop
    End if
    
    !// 边界条件
    Do i = 0, numT
      t = y0 + dble(i) * k
      Lside(i) = exp(-2.d0*t)
      Rside(i) = exp(-1.d0-2.d0*t)
    End do
    u(0,:) = Lside(0:numT)
    u(numH,:) = Rside(0:numT)
    
    !// 初始位移f(x,0) = exp(-x)
    u(:,0) = exp( -1.d0*(xL + [0:numH]*h) )
    
    !// 利用初始位移计算下一个时刻的位移
    Do i = 1, numH-1
      x = xL + dble(i) * h
      u(i,1) = ( 1.d0 - sigma**2 ) * u(i,0) + k*g( x ) + sigma*sigma/2.d0 * ( u(i-1,0) + u(i+1,0) )
    End do
    
    Do j = 2, numT
      Do i = 1, numH-1
        u(i,j) = 2*( 1.d0 - sigma**2 ) * u(i,j-1) + sigma*sigma*( u(i-1,j-1) + u(i+1,j-1) ) - u(i,j-2)
      End do
    End do
    
    !// output
    open( 100, file = 'output.dat' )
    Do i = 0, numT
      t = y0 + dble(i) * k
      Do j = 0, numH
        x = xL + dble(j) * h
        write( 100,'(5f15.8)' ) t, x, u(j,i), exp( -x - 2.d0*t ), abs( exp( -x - 2.d0*t ) - u(j,i) ) / exp( -x - 2.d0*t ) * 100
      End do
    End do
    close( 100 )
    
  End subroutine solve
    
  Real(kind=8) function g( x )  !// 初始速度
    Implicit none
    Real(kind=8), intent(in) :: x
    g = -2.d0*exp( -x )
  End function g
  
End module FwdDeff

Program main
  use FwdDeff
  Implicit none
  call solve
End program main  