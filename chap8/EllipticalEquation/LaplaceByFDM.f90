Module LaplaceByFDM
  Implicit none
  Integer, parameter :: m = 11, n = 11  !// x, y方向的点数
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0
  Real(kind=8), parameter :: ya = 1.d0, yb = 2.d0
  Real(kind=8), parameter :: h = ( xR-xL )/(m-1), k = ( yb-ya )/(n-1)
Contains
  Subroutine solve()
    Implicit none
    Real(kind=8) :: a(m*n,m*n), b(m*n), u(m*n), x(m), y(n), uu(m,n)
    Integer :: i, j
    
    Do i = 1, m
      x(i) = xL + dble(i-1) * h
    End do
    Do i = 1, n
      y(i) = ya + dble(i-1) * k
    End do
    
    !// 求取内部点系数
    Do i = 2, m-1
      Do j = 2, n-1
        a( i+(j-1)*m, (i-1)+(j-1)*m ) = 1.d0 / h**2
        a( i+(j-1)*m, (i+1)+(j-1)*m ) = 1.d0 / h**2
        a( i+(j-1)*m, i+(j-1)*m ) = -2.d0 / h**2 - 2.d0 / k**2
        a( i+(j-1)*m, i+(j-2)*m ) = 1.d0 / k**2
        a( i+(j-1)*m, i+j*m ) = 1.d0 / k**2
        b(i+(j-1)*m) = f( x(i),y(j) )
      End do
    End do
    
    !// 底部与顶部
    Do i = 1, m
      j = 1
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g1( x(i) )
      j = n
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g2( x(i) )
    End do
    
    !// 左侧与右侧
    Do j = 2, n-1
      i = 1
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g3( y(j) )
      i = m
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g4( y(j) )
    End do
    
    u = 0.d0
    call gauss( a, b, u, m*n )
    uu = reshape( u,[m,n] )
    
    !// output
    open( 100, file = 'output.dat' )
    Do i = 1, m
      Do j = 1, n
        write( 100,'(*(f12.4))' ) x(i), y(j), uu(i,j), log( x(i)*x(i) + y(j)*y(j) )
      End do
    End do
    close( 100 )
    
  End subroutine solve
  
  Subroutine gauss( a, b, u, n )
    Implicit none
    Integer, intent(in) :: n
    Real(kind=8), intent(inout) :: a(n,n), b(n), u(n)
    Integer :: i, j, k
    Real(kind=8) :: mult
    
    !// 解方程
    Do j = 1, n - 1
      Do i = j + 1, n
        mult = a(i,j) / a(j,j)
        Do k = j, n 
          a(i,k) = a(i,k) - mult * a(j,k)
        End do 
        b(i) = b(i) - mult * b(j)
      End do 
    End do 
    
    !// 回代
    u = 0.d0
    Do i = n, 1, -1
      Do j = i + 1, n 
        b(i) = b(i) - a(i,j) * u(j)
      End do 
      u(i) = b(i) / a(i,i)
    End do
  End subroutine gauss
  
  !// 右端项以及边界条件
  Real(kind=8) function f( x, y )
    Implicit none
    Real(kind=8), intent(in) :: x, y
    f = 0.d0
  End function f
  
  Real(kind=8) function g1( x )
    Implicit none
    Real(kind=8), intent(in) :: x
    g1 = log( x*x + 1.d0 )
  End function g1
  
  Real(kind=8) function g2( x )
    Implicit none
    Real(kind=8), intent(in) :: x
    g2 = log( x*x + 4.d0 )
  End function g2
  
  Real(kind=8) function g3( y )
    Implicit none
    Real(kind=8), intent(in) :: y
    g3 = 2.d0 * log( y )
  End function g3
  
  Real(kind=8) function g4( y )
    Implicit none
    Real(kind=8), intent(in) :: y
    g4 = log( y*y + 1.d0 )
  End function g4
End module LaplaceByFDM
  
Program main
  use LaplaceByFDM
  Implicit none
  call solve
End program main
!// 原方程: u(x,y) = ln( x^2+y^2 )
!// 偏微分方程
!// △u = 0
!// u(x,0) = ln( x^2+1 )
!// u(x,1) = ln( x^2+4 )
!// u(0,y) = 2*ln(y)
!// u(1,y) = ln(y^2+1)