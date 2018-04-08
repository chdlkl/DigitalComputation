Module LaplaceByFDM
  Implicit none
  Integer, parameter :: m = 11, n = 11  !// x, y����ĵ���
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0
  Real(kind=8), parameter :: ya = 0.d0, yb = 1.d0
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
    
    !// ��ȡ�ڲ���ϵ��
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
    
    !// �ײ��붥��
    Do i = 1, m
      j = 1
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g1( x(i) )
      j = n
      a( i+(j-1)*m, i+(j-1)*m ) = 1.d0
      b( i+(j-1)*m ) = g2( x(i) )
    End do
    
    !// ������Ҳ�
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
        write( 100,'(*(f12.4))' ) x(i), y(j), uu(i,j), exp( -1.d0*x(i)*y(j) ), abs( uu(i,j) - exp( -1.d0*x(i)*y(j) ) ) / exp( -1.d0*x(i)*y(j) ) * 100
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
    
    !// �ⷽ��
    Do j = 1, n - 1
      Do i = j + 1, n
        mult = a(i,j) / a(j,j)
        Do k = j, n 
          a(i,k) = a(i,k) - mult * a(j,k)
        End do 
        b(i) = b(i) - mult * b(j)
      End do 
    End do 
    
    !// �ش�
    u = 0.d0
    Do i = n, 1, -1
      Do j = i + 1, n 
        b(i) = b(i) - a(i,j) * u(j)
      End do 
      u(i) = b(i) / a(i,i)
    End do
  End subroutine gauss
  
  !// �Ҷ����Լ��߽�����
  Real(kind=8) function f( x, y )
    Implicit none
    Real(kind=8), intent(in) :: x, y
    f = exp( -1.d0*x*y ) * ( x*x + y*y )
  End function f
  
  Real(kind=8) function g1( x )  !// ����
    Implicit none
    Real(kind=8), intent(in) :: x
    g1 = 1.d0
  End function g1
  
  Real(kind=8) function g2( x )  !// �ײ�
    Implicit none
    Real(kind=8), intent(in) :: x
    g2 = exp( -1.d0*x )
  End function g2
  
  Real(kind=8) function g3( y )  !// ���
    Implicit none
    Real(kind=8), intent(in) :: y
    g3 = 1.d0
  End function g3
  
  Real(kind=8) function g4( y )  !// �Ҳ�
    Implicit none
    Real(kind=8), intent(in) :: y
    g4 = exp( -1.d0*y )
  End function g4
End module LaplaceByFDM
  
Program main
  use LaplaceByFDM
  Implicit none
  call solve
End program main
!// ԭ����: u(x,y) = exp(-x*y)
!// ƫ΢�ַ���
!// ��u = exp( x*y ) * ( x^2 + y^2 )
!// u(x,0) = 1
!// u(x,1) = exp(-x)
!// u(0,y) = 1
!// u(1,y) = exp(-y)