!// ���ȸ����Գ���������ĸ���
!// ���A'= A����n*n����A�ǶԳƾ������������������x/=0��x'Ax > 0����ƾ���A����������
!// �Գ����������Ƿ������
!// �������A������ʽΪ�����Ӧ������ֵ�ĳ˻�
!// ���ڶԳ��������󣬿�ʹ�ó���˹���ֽⷽ�����ɷֽ�Ϊ: A = R'R������R��һ�������Ǿ���
!// ʹ�ó���˹���ֽ⣬���ڶԳ��������󣬺�һ��ľ�����ȣ�����ֻ��һ�������Ķ���Ԫ�أ�������һ��ļ������ʵ�֣����ҽ���ʹ��һ����ڴ�

Module mod
  Implicit none
  Integer, parameter :: n = 3
  Real(kind=8) :: A(n,n) = [ 4.d0, -2.d0, 2.d0, -2.d0, 2.d0, -4.d0, 2.d0, -4.d0, 11.d0 ]
  Real(kind=8) :: b(n,1) = reshape([6.d0,-10.d0,27.d0],[n,1])
  Real(kind=8) :: R(n,n) = 0.d0, RR(n,n) = 0.d0
  Real(kind=8) :: x(n,1) = 0.d0
Contains
Subroutine GetR ()
  Implicit none
  Integer :: i, j
  Real(kind=8), allocatable :: u(:,:)
  
  u = 0.d0
  Do i = 1, n
    If ( A(i,i) < 0.d0 ) then
      Write ( *,'(1x,g0)' ) '����A���ǶԳ��������󣬳������!'
      stop
    End if
    R(i,i) = sqrt(A(i,i))
    j = n - i
    If ( j >= 1 ) then
      Allocate( u(j,1) )
      u(:,1) = A(i,i+1:n) / R(i,i)
      R(i,i+1:n) = u(:,1)
      A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - matmul( u,transpose(u) )
      Deallocate( u )
    End if
  End do

  RR = transpose(R)

End subroutine GetR

!// Ax = b, A = R'R
!// R'Rx = b, ��R'c = b ���c
!// �����Rx = c���x
Subroutine GetRoot ()
  Implicit none
  Integer :: i, j
  Real(kind=8) :: c(n,1) = 0.d0
  
  !// R'c = b
  Do i = 1, n
    Do j = 1, i - 1
      b(i,1) = b(i,1) - RR(i,j) * c(j,1)
    End do 
    c(i,1) = b(i,1) / RR(i,i)
  End do

  !// ��x: Rx = c
  Do i = n, 1, -1
    Do j = i + 1, n
      c(i,1) = c(i,1) - R(i,j) * x(j,1)
    End do
    x(i,1) = c(i,1) / R(i,i)
  End do

  Write ( *,'(1x,a)' ) 'ԭ���̵Ľ�Ϊ:'
  Do i = 1, n
    Write ( *,'(f12.5)' ) x(i,1)
  End do
End subroutine GetRoot  
End module

Program CholeskyDecomposition
  use mod
  Implicit none
  call GetR ()
  call GetRoot ()
End program CholeskyDecomposition

