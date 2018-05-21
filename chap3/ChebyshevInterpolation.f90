!// Chebyshev interpolation
Module interpolation  
Contains
Subroutine newtdd( x,y,c )  !// 计算牛顿插值系数
  Implicit none 
  Real(kind=8), intent( in ) :: x(:), y(:)
  Real(kind=8), intent( inout ) :: c(:,:)  !// 数组c存储系数矩阵
  Real(kind=8), allocatable :: v(:,:)
  Integer :: i, j, n 

  n = size(x)
  Allocate( v(n,n) )
  v = 0.d0
  Do i = 1, size(x)
    v(i,1) = y(i)
  End Do

  Do j = 2, size(x)
    Do i = 1, size(x) + 1 - j
      v(i,j) = ( v(i+1,j-1) - v(i,j-1) ) / ( x(i+j-1) - x(i) )
    End Do
  End Do 

  c(1,:) = v(1,:)  !// 获取系数
End subroutine newtdd

Subroutine CalInterpolation( x,xx,c,res )
  Implicit none 
  Real(kind=8), intent( in ) :: x(:), xx(:), c(:,:)
  Real(kind=8), intent( inout ) :: res(:,:)
  Real(kind=8), allocatable :: d(:,:)
  Integer :: i, j, m

  m = size(x)
  Allocate( d(m,m) )

  d(1,:) = 1.d0
  Do j = 1, m 
    Do i = 2, m 
      d(i,j) = d(i-1,j) * ( xx(j) - x(i-1) )
    End Do
  End Do
  res = matmul( c,d )
End subroutine CalInterpolation
End module interpolation

Program ChebyshevInterpolation
  use interpolation
  Implicit none 
  Integer, parameter :: n = 50
  Real, parameter :: pi = acos(-1.d0)
  Real(kind=8) :: a = -1.d0, b = 1.d0
  Real(kind=8) :: x(n), xx(n), y(n), c(1,n), res(1,n)
  Integer :: i 

  Do i = 1, n 
    x(i) = (b + a) / 2.d0 + ( b - a ) / 2.d0 * cos( ( 2.d0 * real(i,8) - 1.d0 ) * pi / 2.d0 / real(n,8) )  !// 构造切比雪夫点作为插值基点
    xx(i) = a / 2.d0 + real(i,8) * ( b - a ) / 2.d0 / real(n,8)  !// 待求插值节点
    y(i) = 1.d0 / ( 1.d0 + 12.d0 * x(i) * x(i) )  !// y = 1 / ( 1+12*x^2 )
  End Do 
  Call newtdd( x,y,c )  !// 求解系数c
  Call CalInterpolation( x,xx,c,res )  !// 求解插值节点处的节点值
  Open ( 101, file = '50-PointChebyshev.dat' )
  Do i = 1, n 
    Write ( 101,'(4(2x,g0))' ) x(i), y(i), xx(i), res(1,i)
  End Do 
  Close ( 101 )
End Program ChebyshevInterpolation