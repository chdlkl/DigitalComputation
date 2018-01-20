!// ˵�������ˣ��ڶ����йط������Ĳ��־ʹ˽������½�������ʽ����ڶ��µĵ������֣����ڼ�������
!// -----------ԭ����-----------
!//      2x + 1y + 5z = 5
!//      4x + 4y - 4z = 0
!//      1x + 3y + 1z = 6
!// ---------------------------
!// PA = LU���벿����Ԫ���������ڣ����һ���û�����P��������b
!// ���ڲ�����Ԫ�����A��b��ͬʱ�����б任�ģ�PA = LU �У�A�����б任��b���б任��P��ʵ��
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
  !---------���--------
  !-------Lc = Pb-------
  !--------Ux = c-------
  !---------���x-------
Contains
Subroutine Elimination ( )  !// ��˹��ȥ
  Implicit none 
  Integer :: i, j, k
  Real(kind=8) :: mult, arr(m), arrP(m), P(m,m)

  Write ( *,'(1x,a)' ) '�����в���ǰ����aΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  !// ����P����
  forall ( i = 1:m, j = 1:m, i == j ) P(i,j) = 1.d0
  
  Do j = 1, m - 1
    !//-------------------����Ԫ--------------------
    arr(:) = a(j,:); arrP(:) = P(j,:)
    k = maxloc( abs( a(j:m,j) ), dim = 1 )
    k = k + j - 1
    a(j,:) = a(k,:); P(j,:) = P(k,:)
    a(k,:) = arr(:); P(k,:) = arrP(:)
    !//--------------------------------------------
    !// �����б任��ľ���a
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

  Write ( *,'(1x,a)' ) '�����в��������aΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  Write ( *,'(1x,a)' ) '�����в��������PΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( P(i,j), j = 1, m )
  End do 
  
  !// �Ӿ����б任��ľ���a����ȡ����L
  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
  forall ( i = 1:m, j = 1:m,  i > j ) L(i,j) = a(i,j)
  Write ( *,'(1x,a)' ) '�����в��������LΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
  End do 
  
  !// �Ӿ����б任��ľ���a����ȡ����U
  forall ( i = 1:m, j = 1:m,  i <= j ) U(i,j) = a(i,j)
  Write ( *,'(1x,a)' ) '�����в��������LΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
  End do 
  
  !// ����Pb
  Pb = matmul( P,b )
  
End subroutine Elimination 

Subroutine BackSubstitution ( )  !// �ش�
  Implicit none 
  Integer :: i, j 
  Real(kind=8) :: c(m,1) = 0.d0, x(m,1) = 0.d0
  
  !// ��c: Lc = b
  Do i = 1, m
    Do j = 1, i-1
      Pb(i,1) = Pb(i,1) - L(i,j) * c(j,1)
    End do 
    c(i,1) = Pb(i,1) / L(i,i)
  End do
  
  !// ��x: Ux = c
  Do i = m, 1, -1
    Do j = i + 1, m 
      c(i,1) = c(i,1) - U(i,j) * x(j,1)
    End do 
    x(i,1) = c(i,1) / U(i,i)
  End do 
  
  Write ( *,'(1x,a)' ) 'ԭ���̽�Ϊ��'
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
  
  
!// ������������
!Module mod 
!  Implicit none 
!  Integer, parameter :: m = 4
!  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
!  Real(kind=8) :: b(m,1) = [ 10.d0, 15.d0, 3.d0, 5.d0 ], Pb(m,1) = 0.d0
!  Real(kind=8) :: L(m,m) = 0.d0, U(m,m) = 0.d0
!Contains
!Subroutine Elimination ( )  !// ��˹��ȥ
!  Implicit none 
!  Integer :: i, j, k
!  Real(kind=8) :: mult, arr(m), arrP(m), P(m,m)
!
!  Write ( *,'(1x,a)' ) '�����в���ǰ����aΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  !// ����P����
!  forall ( i = 1:m, j = 1:m, i == j ) P(i,j) = 1.d0
!  
!  Do j = 1, m - 1
!    !//-------------------����Ԫ--------------------
!    arr(:) = a(j,:); arrP(:) = P(j,:)
!    k = maxloc( abs( a(j:m,j) ), dim = 1 )
!    k = k + j - 1
!    a(j,:) = a(k,:); P(j,:) = P(k,:)
!    a(k,:) = arr(:); P(k,:) = arrP(:)
!    !//--------------------------------------------
!    !// �����б任��ľ���a
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
!  Write ( *,'(1x,a)' ) '�����в��������aΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  Write ( *,'(1x,a)' ) '�����в��������PΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( P(i,j), j = 1, m )
!  End do 
!  
!  !// �Ӿ����б任��ľ���a����ȡ����L
!  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
!  forall ( i = 1:m, j = 1:m,  i > j ) L(i,j) = a(i,j)
!  Write ( *,'(1x,a)' ) '�����в��������LΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
!  End do 
!  
!  !// �Ӿ����б任��ľ���a����ȡ����U
!  forall ( i = 1:m, j = 1:m,  i <= j ) U(i,j) = a(i,j)
!  Write ( *,'(1x,a)' ) '�����в��������UΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
!  End do 
!  
!  !// ����Pb
!  Pb = matmul( P,b )
!  
!End subroutine Elimination 
!
!Subroutine BackSubstitution ( )  !// �ش�
!  Implicit none 
!  Integer :: i, j 
!  Real(kind=8) :: c(m,1) = 0.d0, x(m,1) = 0.d0
!  
!  !// ��c: Lc = b
!  Do i = 1, m
!    Do j = 1, i-1
!      Pb(i,1) = Pb(i,1) - L(i,j) * c(j,1)
!    End do 
!    c(i,1) = Pb(i,1) / L(i,i)
!  End do
!  
!  !// ��x: Ux = c
!  Do i = m, 1, -1
!    Do j = i + 1, m 
!      c(i,1) = c(i,1) - U(i,j) * x(j,1)
!    End do 
!    x(i,1) = c(i,1) / U(i,i)
!  End do 
!  
!  Write ( *,'(1x,a)' ) 'ԭ���̽�Ϊ��'
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