!// 说明：至此，第二章有关方程求解的部分就此结束。下节我们正式进入第二章的迭代部分，下期见！！！
!// -----------原方程-----------
!//      2x + 1y + 5z = 5
!//      4x + 4y - 4z = 0
!//      1x + 3y + 1z = 6
!// ---------------------------
!// PA = LU，与部分主元的区别在于，求得一个置换矩阵P，作用于b
!// 而在部分主元里，矩阵A和b是同时进行行变换的，PA = LU 中，A进行行变换，b的行变换由P来实现
Module mod 
  Implicit none 
  Integer, parameter :: m = 3
  Real(kind=8) :: a(m,m) = [ 2.d0, 4.d0, 1.d0, 1.d0, 4.d0, 3.d0, 5.d0, -4.d0, 1.d0 ]
  Real(kind=8) :: b(m,1) = [ 5.d0, 0.d0, 6.d0 ], Pb(m,1) = 0.d0
  Real(kind=8) :: L(m,m) = 0.d0, U(m,m) = 0.d0
  !--------Ax = b-------
  !-------PAx = Pb------
  !--------PA = LU------
  !-------LUx = Pb------
  !---------求解--------
  !-------Lc = Pb-------
  !--------Ux = c-------
  !---------解得x-------
Contains
Subroutine Elimination ( )  !// 高斯消去
  Implicit none 
  Integer :: i, j, k
  Real(kind=8) :: mult, arr(m), arrP(m), P(m,m)

  Write ( *,'(1x,a)' ) '经过行操作前矩阵a为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  !// 构造P矩阵
  forall ( i = 1:m, j = 1:m, i == j ) P(i,j) = 1.d0
  
  Do j = 1, m - 1
    !//-------------------换主元--------------------
    arr(:) = a(j,:); arrP(:) = P(j,:)
    k = maxloc( abs( a(j:m,j) ), dim = 1 )
    k = k + j - 1
    a(j,:) = a(k,:); P(j,:) = P(k,:)
    a(k,:) = arr(:); P(k,:) = arrP(:)
    !//--------------------------------------------
    !// 计算行变换后的矩阵a
    Do i = j + 1, m
      mult = a(i,j) / a(j,j)
      Do k = j, m 
        If ( k == j ) Then 
          a(i,j) = mult
        Else
          a(i,k) = a(i,k) - mult * a(j,k)
        End if 
      End do 
    End do 
  End do 

  Write ( *,'(1x,a)' ) '经过行操作后矩阵a为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  Write ( *,'(1x,a)' ) '经过行操作后矩阵P为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( P(i,j), j = 1, m )
  End do 
  
  !// 从经过行变换后的矩阵a中提取矩阵L
  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
  forall ( i = 1:m, j = 1:m,  i > j ) L(i,j) = a(i,j)
  Write ( *,'(1x,a)' ) '经过行操作后矩阵L为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
  End do 
  
  !// 从经过行变换后的矩阵a中提取矩阵U
  forall ( i = 1:m, j = 1:m,  i <= j ) U(i,j) = a(i,j)
  Write ( *,'(1x,a)' ) '经过行操作后矩阵L为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
  End do 
  
  !// 计算Pb
  Pb = matmul( P,b )
  
End subroutine Elimination 

Subroutine BackSubstitution ( )  !// 回代
  Implicit none 
  Integer :: i, j 
  Real(kind=8) :: c(m,1) = 0.d0, x(m,1) = 0.d0
  
  !// 求c: Lc = b
  Do i = 1, m
    Do j = 1, i-1
      Pb(i,1) = Pb(i,1) - L(i,j) * c(j,1)
    End do 
    c(i,1) = Pb(i,1) / L(i,i)
  End do
  
  !// 求x: Ux = c
  Do i = m, 1, -1
    Do j = i + 1, m 
      c(i,1) = c(i,1) - U(i,j) * x(j,1)
    End do 
    x(i,1) = c(i,1) / U(i,i)
  End do 
  
  Write ( *,'(1x,a)' ) '原方程解为：'
  Do i = 1, m 
    Write ( *,'(f12.5)' ) x(i,1) 
  End do 
  
End subroutine BackSubstitution
End module mod 


Program GaussianElimination
  Use mod 
  Implicit none 
  call Elimination ( )
  call BackSubstitution ( )
End program GaussianElimination
  
  
!// 测试其他方程
!Module mod 
!  Implicit none 
!  Integer, parameter :: m = 4
!  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
!  Real(kind=8) :: b(m,1) = [ 10.d0, 15.d0, 3.d0, 5.d0 ], Pb(m,1) = 0.d0
!  Real(kind=8) :: L(m,m) = 0.d0, U(m,m) = 0.d0
!Contains
!Subroutine Elimination ( )  !// 高斯消去
!  Implicit none 
!  Integer :: i, j, k
!  Real(kind=8) :: mult, arr(m), arrP(m), P(m,m)
!
!  Write ( *,'(1x,a)' ) '经过行操作前矩阵a为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  !// 构造P矩阵
!  forall ( i = 1:m, j = 1:m, i == j ) P(i,j) = 1.d0
!  
!  Do j = 1, m - 1
!    !//-------------------换主元--------------------
!    arr(:) = a(j,:); arrP(:) = P(j,:)
!    k = maxloc( abs( a(j:m,j) ), dim = 1 )
!    k = k + j - 1
!    a(j,:) = a(k,:); P(j,:) = P(k,:)
!    a(k,:) = arr(:); P(k,:) = arrP(:)
!    !//--------------------------------------------
!    !// 计算行变换后的矩阵a
!    Do i = j + 1, m
!      mult = a(i,j) / a(j,j)
!      Do k = j, m 
!        If ( k == j ) Then 
!          a(i,j) = mult
!        Else
!          a(i,k) = a(i,k) - mult * a(j,k)
!        End if 
!      End do 
!    End do 
!  End do 
!
!  Write ( *,'(1x,a)' ) '经过行操作后矩阵a为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  Write ( *,'(1x,a)' ) '经过行操作后矩阵P为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( P(i,j), j = 1, m )
!  End do 
!  
!  !// 从经过行变换后的矩阵a中提取矩阵L
!  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
!  forall ( i = 1:m, j = 1:m,  i > j ) L(i,j) = a(i,j)
!  Write ( *,'(1x,a)' ) '经过行操作后矩阵L为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
!  End do 
!  
!  !// 从经过行变换后的矩阵a中提取矩阵U
!  forall ( i = 1:m, j = 1:m,  i <= j ) U(i,j) = a(i,j)
!  Write ( *,'(1x,a)' ) '经过行操作后矩阵U为：'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
!  End do 
!  
!  !// 计算Pb
!  Pb = matmul( P,b )
!  
!End subroutine Elimination 
!
!Subroutine BackSubstitution ( )  !// 回代
!  Implicit none 
!  Integer :: i, j 
!  Real(kind=8) :: c(m,1) = 0.d0, x(m,1) = 0.d0
!  
!  !// 求c: Lc = b
!  Do i = 1, m
!    Do j = 1, i-1
!      Pb(i,1) = Pb(i,1) - L(i,j) * c(j,1)
!    End do 
!    c(i,1) = Pb(i,1) / L(i,i)
!  End do
!  
!  !// 求x: Ux = c
!  Do i = m, 1, -1
!    Do j = i + 1, m 
!      c(i,1) = c(i,1) - U(i,j) * x(j,1)
!    End do 
!    x(i,1) = c(i,1) / U(i,i)
!  End do 
!  
!  Write ( *,'(1x,a)' ) '原方程解为：'
!  Do i = 1, m 
!    Write ( *,'(f12.5)' ) x(i,1) 
!  End do 
!  
!End subroutine BackSubstitution
!End module mod 
!
!
!Program GaussianElimination
!  Use mod 
!  Implicit none 
!  call Elimination ( )
!  call BackSubstitution ( )
!End program GaussianElimination