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
Subroutine Jacobi ()  !// Jacobi�������
  Implicit none
  Real(kind=8) :: A(n,n), InvD(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = ��ʼ����
  !// Xk+1 = InvD[ b - (L+U)Xk ], k = 0, 1, 2,...
  !// InvDΪϵ������Խ���Ԫ�ص������
  !// bΪ�Ҷ���
  !// LΪϵ������������ǲ��֣�ע����LU�ֽ��е�L��ͬ
  !// UΪϵ������������ǲ��֣�ע����LU�ֽ��е�U��ͬ
  !// XkΪǰһ�μ�������Ľ��
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
  
  x0 = 0.d0  !// ��ʼ������
  
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

Subroutine Gauss_Seidel ()  !// ��˹-���¶��������
  Implicit none 
  Real(kind=8) :: A(n,n), InvD(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = ��ʼ����
  !// Xk+1 = InvD[ b - U*Xk - L*Xk+1 ], k = 0, 1, 2,...
  !//--------------------------------------------------
  Open ( newunit = fileid, file = 'IterationData.txt' )
  Read ( fileid, * ) A  !// ����ע��һ�£������뽫�����е�ϡ�������Ҷ���д���ļ��н��ж�ȡ�����߿ɸ����Լ�������ʵ��޸�
  Read ( fileid, * ) b
  Close ( fileid )
  !//--------------------------------------------------
  
  InvD = 0.d0
  forall ( i = 1:n, j = 1:n, i == j ) InvD(i,j) = 1.d0 / A(i,j)
  
  L = 0.d0 
  forall ( i = 1:n, j = 1:n, i > j ) L(i,j) = A(i,j)
  
  U = 0.d0
  forall ( i = 1:n, j = 1:n, i < j ) U(i,j) = A(i,j)
  
  x0 = 0.d0  !// ��ʼ������
  
  i = 1
  Do 
    !//------------------------------------------
    tmp = b - matmul( (L+U),x0 )  !// X0Ϊʽ�е�Xk
    x = matmul( InvD, tmp )  !// �ȼ�����Ҳ��Xk+1���Ҳ��Xk+1��Jacobi����õ�
    !//------------------------------------------
    tmp = b - matmul( U, x0 ) - matmul( L, x )
    x = matmul( InvD, tmp )  !// �����������Xk+1
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

Subroutine SOR ()  !// �������ɳڵ������
  Implicit none 
  Real(kind=8) :: A(n,n), D(n,n), L(n,n), U(n,n)
  Real(kind=8) :: b(n,1), x0(n,1), x(n,1), tmp(n,1)
  Real(kind=8) :: LD(n,n), InvLD(n,n)  !// LD = wL + D
  Real(kind=8), parameter :: w = 1.1d0  !// wΪ�ɳ����ӣ�w����0ʱ���ӿ����������ɳڣ���С��0ʱ����������
  Integer :: fileid, i, j
  !// Jacobi
  !// X0 = ��ʼ����
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
  
  x0 = 0.d0  !// ��ʼ������
  !// ����wL + D����
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

Subroutine Inv ( aa, b, n )  !// �������
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
!// ���������ֵ��������������ϸ�Խ�ռ�ž��󣬶�����������|Aij|>sum(|Aij|,i/=j)
!// ���n*n����Aʱ�ϸ�ĶԽ�ռ�ž�����(1):A�Ƿ��������
!// (2)����������b�ͳ�ʼ���ƣ���Ax=bӦ������ĵ�����������������(Ψһ)�⡣
!// ע�⣺�ϸ�Խ�ռ�Ž�����һ�����������������Խ�ռ��ʱ����Ȼ����������
!// �������������������ϸ�Խ�ռ�еģ����������һ���󣬿���ʹ��PAeqLU2.0.f90�����û�����P�Ĵ���Σ�����һ����ת���ɶԽ�����Ԫ�����ľ���������