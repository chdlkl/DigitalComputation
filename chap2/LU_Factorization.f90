!//---------------原方程-------------
!//     x + 2y - z = 3
!//    2x + y -2z = 3
!//    -3x + y + z = -6
!//----------------------------------
!// 本代实现了LU分解，但是并不是所有的矩阵都能进行LU分解
!// 后面会介绍PA=LU，可以解决所有问题
!// 为什么要进行LU分解
!// 使用LU分解方法代替高斯消去法的主要原因是由于如下系统的存在
!// Ax=B1, Ax=B2, Ax=B3,...,Ax=Bk。一般称A为结构矩阵
!// 经典高斯消去法对每一个Ax=B都要分解，但是因为系统有相同的A，LU分解将Bk独立出来，只进行一次LU分解即可
Module mod 
  Implicit none 
  Integer, parameter :: m = 3  !// 方程的阶数
  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, -3.d0, 2.d0, 1.d0, 1.d0, -1.d0, -2.d0, 1.d0 ]
  Real(kind=8) :: origin_b(m) = [ 3.d0, 3.d0, -6.d0 ] 
  Real(kind=8) :: c(m) = 0.d0, x(m) = 0.d0 
  Real(kind=8) :: U(m,m) = 0.d0, L(m,m) = 0.d0
  !// origin_b存放初始右端项，用来计算Lc = b,计算c
  !// 然后利用Ux = c，求解x
Contains
Subroutine GetLU ( )  
  Implicit none 
  Integer :: i, j, k 
  Real(kind=8), parameter :: eps = 1.d-5  !// 当主元小于这个数时，程序退出
  Real(kind=8) :: mult

  Write ( *,'(1x,a)' ) '原系数矩阵a为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  !// 对下三角矩阵对角线赋值
  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
  
  !// 计算上三角矩阵和下三角矩阵
  U = a
  Do j = 1, m - 1
    If ( abs(U(j,j)) < eps ) Then
      Write ( *,'(1x,a)' ) ' The pivot is zero!'
      stop 
    End if 
    Do i = j + 1, m
      mult = U(i,j) / U(j,j)
      L(i,j) = mult
      Do k = j, m 
        U(i,k) = U(i,k) - mult * U(j,k)
      End do 
    End do 
  End do 
  
  !// 得到上三角矩阵
  Write ( *,'(1x,a)' ) '上三角矩阵为：'
  Do i = 1, m
    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
  End do
  
  Write ( *,'(1x,a)' ) '下三角矩阵为：'
  Do i = 1, m
    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
  End do
  
End subroutine GetLU 

Subroutine BackSubstitution ( )
  Implicit none 
  Integer :: i, j 
  
  !// 求c: Lc = origin_b
  Do i = 1, m
    Do j = 1, i-1
      origin_b(i) = origin_b(i) - L(i,j) * c(j)
    End do 
    c(i) = origin_b(i) / L(i,i)
  End do
  
  !// 求x: Ux = c
  Do i = m, 1, -1
    Do j = i + 1, m 
      c(i) = c(i) - U(i,j) * x(j)
    End do 
    x(i) = c(i) / U(i,i)
  End do 
  
  Write ( *,'(1x,a)' ) '原方程解为：'
  Do i = 1, m 
    Write ( *,'(f12.5)' ) x(i) 
  End do 
  
End subroutine BackSubstitution

End module mod 


Program LU_Factorization
  Use mod 
  Implicit none 
  call GetLU ( )
  call BackSubstitution ( )
End program LU_Factorization

!// 对另一方程组进行验证。x=[1,2,3,4]  
!//---------------原方程-------------
!//     x + y + z + w = 10
!//    2x + 3y + z + w = 15
!//    3x - y + 2z - w = 3
!//    4x + y -3z + 2w = 5
!//----------------------------------
!Module mod 
!  Implicit none 
!  Integer, parameter :: m = 4  !// 方程的阶数
!  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
!  Real(kind=8) :: origin_b(m) = [ 10.d0, 15.d0, 3.d0, 5.d0 ] 
!  Real(kind=8) :: c(m) = 0.d0, x(m) = 0.d0 
!  Real(kind=8) :: U(m,m) = 0.d0, L(m,m) = 0.d0
!  !// origin_b存放初始右端项，用来计算Lc = b,计算c
!  !// 然后利用Ux = c，求解x
!Contains
!Subroutine GetLU ( )  
!  Implicit none 
!  Integer :: i, j, k 
!  Real(kind=8), parameter :: eps = 1.d-5  !// 当主元小于这个数时，程序退出
!  Real(kind=8) :: mult
!
!  Write ( *,'(1x,a)' ) '原系数矩阵a为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  !// 对下三角矩阵对角线赋值
!  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
!  
!  !// 计算上三角矩阵和下三角矩阵
!  U = a
!  Do j = 1, m - 1
!    If ( abs(U(j,j)) < eps ) Then
!      Write ( *,'(1x,a)' ) ' The pivot is zero!'
!      stop 
!    End if 
!    Do i = j + 1, m
!      mult = U(i,j) / U(j,j)
!      L(i,j) = mult
!      Do k = j, m 
!        U(i,k) = U(i,k) - mult * U(j,k)
!      End do 
!    End do 
!  End do 
!  
!  !// 得到上三角矩阵
!  Write ( *,'(1x,a)' ) '上三角矩阵为：'
!  Do i = 1, m
!    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
!  End do
!  
!  Write ( *,'(1x,a)' ) '下三角矩阵为：'
!  Do i = 1, m
!    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
!  End do
!  
!End subroutine GetLU 
!
!Subroutine BackSubstitution ( )
!  Implicit none 
!  Integer :: i, j 
!  
!  !// 求c: Lc = origin_b
!  Do i = 1, m
!    Do j = 1, i-1
!      origin_b(i) = origin_b(i) - L(i,j) * c(j)
!    End do 
!    c(i) = origin_b(i) / L(i,i)
!  End do
!  
!  !// 求x: Ux = c
!  Do i = m, 1, -1
!    Do j = i + 1, m 
!      c(i) = c(i) - U(i,j) * x(j)
!    End do 
!    x(i) = c(i) / U(i,i)
!  End do 
!  
!  Write ( *,'(1x,a)' ) '原方程解为：'
!  Do i = 1, m 
!    Write ( *,'(f12.5)' ) x(i) 
!  End do 
!  
!End subroutine BackSubstitution
!
!End module mod 
!
!
!Program LU_Factorization
!  Use mod 
!  Implicit none 
!  call GetLU ( )
!  call BackSubstitution ( )
!End program LU_Factorization
!
