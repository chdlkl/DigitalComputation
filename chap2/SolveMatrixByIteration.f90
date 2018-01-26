!// [ 3 -1 0  0 0 0.5 ] [x1]   [2.5]
!// [ -1 3 -1 0 0.5 0 ] [x2]   [1.5]
!// [ 0 -1 3 -1  0  0 ] [x3] = [1.0]
!// [ 0  0 -1  3 -1 0 ] [x4]   [1.0]
!// [ 0 0.5 0 -1 3 -1 ] [x5]   [1.5]
!// [ 0.5 0 0 0 -1 3  ] [x6]   [2.5]
Module Iteration
  Implicit none
  Real(kind=8) :: eps = 1d-12, err
  Integer, parameter :: maxInteration = 50, n = 6
Contains
Subroutine Jacobi ()  !// Jacobi迭代求解
  Implicit none
  Real(kind=8) :: A(n,n), InvD(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = 初始向量
  !// Xk+1 = InvD[ b - (L+U)Xk ], k = 0, 1, 2,...
  !// InvD为系数矩阵对角线元素的逆矩阵
  !// b为右端项
  !// L为系数矩阵的下三角部分，注意与LU分解中的L不同
  !// U为系数矩阵的上三角部分，注意与LU分解中的U不同
  !// Xk为前一次计算出来的结果
  Open ( newunit = fileid, file = 'IterationData.txt' )
  Read ( fileid, * ) A
  Read ( fileid, * ) b
  Close ( fileid )
  
  InvD = 0.d0
  forall ( i = 1:n, j = 1:n, i == j ) InvD(i,j) = 1.d0 / A(i,j)
  
  L = 0.d0 
  forall ( i = 1:n, j = 1:n, i > j ) L(i,j) = A(i,j)
  
  U = 0.d0
  forall ( i = 1:n, j = 1:n, i < j ) U(i,j) = A(i,j)
  
  x0 = 0.d0  !// 初始化向量
  
  i = 1
  Do 
    tmp = b - matmul( (L+U),x0 )
    x = matmul( InvD, tmp )
    i = i + 1
    err = maxval( abs(x-x0) )
    If ( i > maxInteration .or. err < eps ) exit
    x0 = x
  End do
  Write ( *,'(1x,A)' ) "Jacobi solution: "
  Write ( *,'(*(f9.6))' ) x
  Write ( *,'(1x,A,I3)' ) "The iterations of Jacobi is ", i
  
End subroutine Jacobi

Subroutine Gauss_Seidel ()  !// 高斯-赛德尔迭代求解
  Implicit none 
  Real(kind=8) :: A(n,n), InvD(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = 初始向量
  !// Xk+1 = InvD[ b - U*Xk - L*Xk+1 ], k = 0, 1, 2,...
  !//--------------------------------------------------
  Open ( newunit = fileid, file = 'IterationData.txt' )
  Read ( fileid, * ) A  !// 这里注意一下，本代码将矩阵中的稀疏矩阵和右端项写入文件中进行读取，读者可根据自己的情况适当修改
  Read ( fileid, * ) b
  Close ( fileid )
  !//--------------------------------------------------
  
  InvD = 0.d0
  forall ( i = 1:n, j = 1:n, i == j ) InvD(i,j) = 1.d0 / A(i,j)
  
  L = 0.d0 
  forall ( i = 1:n, j = 1:n, i > j ) L(i,j) = A(i,j)
  
  U = 0.d0
  forall ( i = 1:n, j = 1:n, i < j ) U(i,j) = A(i,j)
  
  x0 = 0.d0  !// 初始化向量
  
  i = 1
  Do 
    !//------------------------------------------
    tmp = b - matmul( (L+U),x0 )  !// X0为式中的Xk
    x = matmul( InvD, tmp )  !// 先计算出右侧的Xk+1，右侧的Xk+1由Jacobi计算得到
    !//------------------------------------------
    tmp = b - matmul( U, x0 ) - matmul( L, x )
    x = matmul( InvD, tmp )  !// 最后计算出左侧的Xk+1
    !//------------------------------------------
    i = i + 1
    err = maxval( abs(x-x0) )
    If ( i > maxInteration .or. err < eps ) exit
    x0 = x
  End do
  Write ( *,'(1x,A)' ) "Gauss_Seidel solution: "
  Write ( *,'(*(f9.6))' ) x
  Write ( *,'(1x,A,I3)' ) "The iterations of Gauss_Seidel is ", i
  
End subroutine Gauss_Seidel

Subroutine SOR ()  !// 连续过松弛迭代求解
  Implicit none 
  Real(kind=8) :: A(n,n), D(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Real(kind=8) :: LD(n,n), InvLD(n,n)  !// LD = wL + D
  Real(kind=8), parameter :: w = 1.1d0  !// w为松弛因子，w大于0时，加快收敛（过松弛），小于0时，减缓收敛
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = 初始向量
  !// Xk+1 = Inv( wL + D ) * [ (1-w)*D*Xk - w*U*Xk ] + w*Inv( w*L + D )*b, k = 0, 1, 2,...
  Open ( newunit = fileid, file = 'IterationData.txt' )
  Read ( fileid, * ) A
  Read ( fileid, * ) b
  Close ( fileid )
  
  D = 0.d0
  forall ( i = 1:n, j = 1:n, i == j ) D(i,j) = A(i,j)
  
  L = 0.d0 
  forall ( i = 1:n, j = 1:n, i > j ) L(i,j) = A(i,j)
  
  U = 0.d0
  forall ( i = 1:n, j = 1:n, i < j ) U(i,j) = A(i,j)
  
  x0 = 0.d0  !// 初始化向量
  !// 计算wL + D的逆
  LD = w * L + D
  InvLD = 0.d0
  call Inv ( LD, InvLD, n )
  
  i = 1
  Do 
    !//------------------------------------------
    tmp = matmul ( ( 1.d0 - w )*D, x0 ) - matmul( w*U, x0 )  !// (1-w)*D*Xk - w*U*Xk
    x = matmul( InvLD, tmp ) + matmul( w*InvLD, b )  !// Xk+1 = Inv( wL + D ) * [ (1-w)*D*Xk - w*U*Xk ] + w*Inv( w*L + D )*b
    !//------------------------------------------
    i = i + 1
    err = maxval( abs(x-x0) )
    If ( i > maxInteration .or. err < eps ) exit
    x0 = x
  End do
  Write ( *,'(1x,A)' ) "SOR solution: "
  Write ( *,'(*(f9.6))' ) x
  Write ( *,'(1x,A,I3)' ) "The iterations of SOR ", i
  
End subroutine SOR

Subroutine Inv ( aa, b, n )  !// 求逆矩阵
  Implicit none
  Integer :: n,i,j,k
  Real(kind=8) :: aa(n,n), b(n,n), a(n,n)
  
  a = aa
  Do i = 1, n
    b(i,i) = 1.d0
  End do
  
  Do i = 1, n
    b(i,:) = b(i,:) / a(i,i)
    a(i,i:n) = a(i,i:n) / a(i,i)
    Do j = i + 1, n
      Do k = 1, n
        b(j,k) = b(j,k) - b(i,k)*a(j,i)
      End do
      a(j,i:n) = a(j,i:n) - a(i,i:n)*a(j,i)
    End do
  End do
  
  Do i = n, 1, -1
    Do j = i - 1, 1, -1
      Do k = 1, n
        b(j,k) = b(j,k) - b(i,k)*a(j,i)
      End do
    End do
  End do
  
End subroutine Inv

End module Iteration
  
  
Program SolveByIteration
  Use Iteration
  Implicit none 
  call Jacobi ()
  call Gauss_Seidel ()
  call SOR ()
End program SolveByIteration
!// 上述的三种迭代方法，对于严格对角占优矩阵，都可以收敛。|Aij|>sum(|Aij|,i/=j)
!// 如果n*n矩阵A时严格的对角占优矩阵，则(1):A是非奇异矩阵
!// (2)对所有向量b和初始估计，对Ax=b应用上面的迭代方法都会收敛到(唯一)解。
!// 注意：严格对角占优仅仅是一个充分条件，不满足对角占优时，依然可能收敛。
!// 本代码所给的例子是严格对角占有的，如果对于任一矩阵，可以使用PAeqLU2.0.f90中求置换矩阵P的代码段，将任一矩阵转化成对角线上元素最大的矩阵进行求解