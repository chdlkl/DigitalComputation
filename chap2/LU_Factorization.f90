!//---------------ԭ����-------------
!//     x + 2y - z = 3
!//    2x + y -2z = 3
!//    -3x + y + z = -6
!//----------------------------------
!// ����ʵ����LU�ֽ⣬���ǲ��������еľ����ܽ���LU�ֽ�
!// ��������PA=LU�����Խ����������
!// ΪʲôҪ����LU�ֽ�
!// ʹ��LU�ֽⷽ�������˹��ȥ������Ҫԭ������������ϵͳ�Ĵ���
!// Ax=B1, Ax=B2, Ax=B3,...,Ax=Bk��һ���AΪ�ṹ����
!// �����˹��ȥ����ÿһ��Ax=B��Ҫ�ֽ⣬������Ϊϵͳ����ͬ��A��LU�ֽ⽫Bk����������ֻ����һ��LU�ֽ⼴��
Module mod 
  Implicit none 
  Integer, parameter :: m = 3  !// ���̵Ľ���
  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, -3.d0, 2.d0, 1.d0, 1.d0, -1.d0, -2.d0, 1.d0 ]
  Real(kind=8) :: origin_b(m) = [ 3.d0, 3.d0, -6.d0 ] 
  Real(kind=8) :: c(m) = 0.d0, x(m) = 0.d0 
  Real(kind=8) :: U(m,m) = 0.d0, L(m,m) = 0.d0
  !// origin_b��ų�ʼ�Ҷ����������Lc = b,����c
  !// Ȼ������Ux = c�����x
Contains
Subroutine GetLU ( )  
  Implicit none 
  Integer :: i, j, k 
  Real(kind=8), parameter :: eps = 1.d-5  !// ����ԪС�������ʱ�������˳�
  Real(kind=8) :: mult

  Write ( *,'(1x,a)' ) 'ԭϵ������aΪ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
  End do 
  
  !// �������Ǿ���Խ��߸�ֵ
  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
  
  !// ���������Ǿ���������Ǿ���
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
  
  !// �õ������Ǿ���
  Write ( *,'(1x,a)' ) '�����Ǿ���Ϊ��'
  Do i = 1, m
    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
  End do
  
  Write ( *,'(1x,a)' ) '�����Ǿ���Ϊ��'
  Do i = 1, m
    Write ( *,'(*(f12.5))' ) ( L(i,j), j = 1, m )
  End do
  
End subroutine GetLU 

Subroutine BackSubstitution ( )
  Implicit none 
  Integer :: i, j 
  
  !// ��c: Lc = origin_b
  Do i = 1, m
    Do j = 1, i-1
      origin_b(i) = origin_b(i) - L(i,j) * c(j)
    End do 
    c(i) = origin_b(i) / L(i,i)
  End do
  
  !// ��x: Ux = c
  Do i = m, 1, -1
    Do j = i + 1, m 
      c(i) = c(i) - U(i,j) * x(j)
    End do 
    x(i) = c(i) / U(i,i)
  End do 
  
  Write ( *,'(1x,a)' ) 'ԭ���̽�Ϊ��'
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

!// ����һ�����������֤��x=[1,2,3,4]  
!//---------------ԭ����-------------
!//     x + y + z + w = 10
!//    2x + 3y + z + w = 15
!//    3x - y + 2z - w = 3
!//    4x + y -3z + 2w = 5
!//----------------------------------
!Module mod 
!  Implicit none 
!  Integer, parameter :: m = 4  !// ���̵Ľ���
!  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
!  Real(kind=8) :: origin_b(m) = [ 10.d0, 15.d0, 3.d0, 5.d0 ] 
!  Real(kind=8) :: c(m) = 0.d0, x(m) = 0.d0 
!  Real(kind=8) :: U(m,m) = 0.d0, L(m,m) = 0.d0
!  !// origin_b��ų�ʼ�Ҷ����������Lc = b,����c
!  !// Ȼ������Ux = c�����x
!Contains
!Subroutine GetLU ( )  
!  Implicit none 
!  Integer :: i, j, k 
!  Real(kind=8), parameter :: eps = 1.d-5  !// ����ԪС�������ʱ�������˳�
!  Real(kind=8) :: mult
!
!  Write ( *,'(1x,a)' ) 'ԭϵ������aΪ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m )
!  End do 
!  
!  !// �������Ǿ���Խ��߸�ֵ
!  forall ( i = 1:m, j = 1:m, i == j ) L(i,j) = 1.d0
!  
!  !// ���������Ǿ���������Ǿ���
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
!  !// �õ������Ǿ���
!  Write ( *,'(1x,a)' ) '�����Ǿ���Ϊ��'
!  Do i = 1, m
!    Write ( *,'(*(f12.5))' ) ( U(i,j), j = 1, m )
!  End do
!  
!  Write ( *,'(1x,a)' ) '�����Ǿ���Ϊ��'
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
!  !// ��c: Lc = origin_b
!  Do i = 1, m
!    Do j = 1, i-1
!      origin_b(i) = origin_b(i) - L(i,j) * c(j)
!    End do 
!    c(i) = origin_b(i) / L(i,i)
!  End do
!  
!  !// ��x: Ux = c
!  Do i = m, 1, -1
!    Do j = i + 1, m 
!      c(i) = c(i) - U(i,j) * x(j)
!    End do 
!    x(i) = c(i) / U(i,i)
!  End do 
!  
!  Write ( *,'(1x,a)' ) 'ԭ���̽�Ϊ��'
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
