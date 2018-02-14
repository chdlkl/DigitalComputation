!// 首先给出对称正定矩阵的概念
!// 如果A'= A，则n*n矩阵A是对称矩阵。如果对于所有向量x/=0，x'Ax > 0，则称矩阵A是正定矩阵
!// 对称正定矩阵是非奇异的
!// 任意矩阵A的行列式为矩阵对应的特征值的乘积
!// 对于对称正定矩阵，可使用楚列斯基分解方法，可分解为: A = R'R，其中R是一个上三角矩阵
!// 使用楚列斯基分解，对于对称正定矩阵，和一般的矩阵相比，它们只有一半数量的独立元素，可以用一半的计算代价实现，并且仅仅使用一半的内存

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
      Write ( *,'(1x,g0)' ) '矩阵A不是对称正定矩阵，程序结束!'
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
!// R'Rx = b, 令R'c = b 求出c
!// 最后用Rx = c求出x
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

  !// 求x: Rx = c
  Do i = n, 1, -1
    Do j = i + 1, n
      c(i,1) = c(i,1) - R(i,j) * x(j,1)
    End do
    x(i,1) = c(i,1) / R(i,i)
  End do

  Write ( *,'(1x,a)' ) '原方程的解为:'
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

