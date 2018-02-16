!// 在病态矩阵问题上，共轭梯度法比部分主元的高斯消去法还要差
!// 这种情况下，可通过预条件得到缓解，主要是将问题转化为良态矩阵系统，再实施CG方法
!// 此代码的预条件因子M= ( D + wL ) * Inv(D) * ( D + wU )，为高斯-赛德尔预条件因子
!// A = L + D + U，L为A的下三角，U为A的上三角，D为A的对角线元素
Program ConjugateGradientMethods
  Implicit none
  Integer :: i, j
  Integer, parameter :: n = 3
  Real(kind=8) :: A(n,n) = [ 4.d0, -2.d0, 2.d0, -2.d0, 2.d0, -4.d0, 2.d0, -4.d0, 11.d0 ]
  Real(kind=8) :: b(n,1) = [ 6.d0, -10.d0, 27.d0 ]
  Real(kind=8) :: alpha = 0.d0, beta = 0.d0
  Real(kind=8) :: x(n,1) = 0.d0, d(n,1) = 0.d0, r(n,1) = 0.d0, z(n,1) = 0.d0
  Real(kind=8) :: rtmp(n,1) = 0.d0, tmp(1,1) = 0.d0, ztmp(n,1) = 0.d0
  Real(kind=8) :: M(n,n) = 0.d0
  
  call GetM ( A, M, n )
  
  r = b - matmul( A,x ) !// r0 = b - Ax0
  z = matmul( M,r )  !// z0 = Inv(M) * r0
  d = z  !// d0 = z0

  Do i = 0, n - 1
    If ( all(r<0.0) ) then
      exit
    End if
    rtmp = r
    ztmp = z
    tmp = matmul( transpose(rtmp),ztmp ) / matmul( matmul( transpose(d),A ), d )
    alpha = tmp(1,1)
    x = x + alpha * d
    r = r - alpha * matmul( A,d )
    z = matmul( M,r )
    tmp = matmul( transpose(r),z ) / matmul( transpose(rtmp),ztmp )
    beta = tmp(1,1)
    d = z + beta * d
  End do
  
  Do i = 1, n
    Write ( *,'(f9.3)' ) x(i,1)
  End do
End program ConjugateGradientMethods
  
Subroutine GetM ( A, M, n )
  Implicit none
  Integer :: i, j
  Integer, intent(in) :: n 
  Real(kind=8), intent(inout) :: A(n,n), M(n,n)
  Real(kind=8) :: L(n,n), D(n,n), U(n,n)
  Real(kind=8), parameter :: w = 1.d0  !// w为高斯-赛德尔预条件因子，w是0-2之间的一个实数
  
  forall ( i = 1:n, j = 1:n, i == j ) D(i,j) = A(i,j)
  forall ( i = 1:n, j = 1:n, i == j ) D(i,j) = 1.d0 / D(i,j)
  forall ( i = 1:n, j = 1:n, i < j ) U(i,j) = A(i,j)
  forall ( i = 1:n, j = 1:n, i > j ) L(i,j) = A(i,j)
  !// M= ( D + wL ) * Inv(D) * ( D + wU )
  M = D + w * L
  M = matmul( M, D )
  M = matmul( M, D+w*U )
End subroutine GetM