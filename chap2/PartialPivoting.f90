!//---------------ԭ����-------------
!//     x + y + z + w = 10
!//    2x + 3y + z + w = 15
!//    3x - y + 2z - w = 3
!//    4x + y -3z + 2w = 5
!//----------------------------------
!// ��������Ը�˹��Ԫ�����иĽ������ò�����Ԫ��˼·

Module mod 
  Implicit none 
  Integer, parameter :: m = 4
  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
  Real(kind=8) :: b(m) = [ 10.d0, 15.d0, 3.d0, 5.d0 ]
  Real(kind=8) :: x(m) = 0.d0 
Contains
Subroutine Elimination ( )  !// ��˹��ȥ
  Implicit none 
  Integer :: i, j, k
  Real(kind=8) :: mult, arr(m), arr_b

  Write ( *,'(1x,a)' ) '������ȥǰ��������Ҷ���Ϊ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m ), b(i)
  End do 
  
  Do j = 1, m - 1
    !//---------------����Ԫ--------------
    arr(:) = a(j,:); arr_b = b(j)
    k = maxloc( abs( a(j:m,j) ), dim = 1 )
    k = k + j - 1
    a(j,:) = a(k,:); b(j) = b(k)
    a(k,:) = arr(:); b(k) = arr_b
    !//----------------------------------
    Do i = j + 1, m
      mult = a(i,j) / a(j,j)
      Do k = j, m 
        a(i,k) = a(i,k) - mult * a(j,k)
      End do 
      b(i) = b(i) - mult * b(j)
    End do 
  End do 

  Write ( *,'(1x,a)' ) '������ȥ����������Ҷ���Ϊ��'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m ), b(i)
  End do 
  
End subroutine Elimination 

Subroutine BackSubstitution ( )
  Implicit none 
  Integer :: i, j 
  
  Do i = m, 1, -1
    Do j = i + 1, m 
      b(i) = b(i) - a(i,j) * x(j)
    End do 
    x(i) = b(i) / a(i,i)
  End do 

  Write ( *,'(1x,a)' ) 'ԭ���̽�Ϊ��'
  Do i = 1, m 
    Write ( *,'(f12.5)' ) x(i) 
  End do 
  
End subroutine BackSubstitution

End module mod 


Program GaussianElimination
  Use mod 
  Implicit none 
  call Elimination ( )
  call BackSubstitution ( )
End program GaussianElimination
 
!// �ڶ����ִ��뾭���˹��Ԫ��ʧЧ�������ǸĽ���Ĵ��루������Ԫ���������Գɹ�
!Module mod 
!  Implicit none 
!  Integer, parameter :: m = 3
!  Real(kind=8) :: a(m,m) = [ 0, 1, 1, 1, 1, 1, 1, 1, -1 ]
!  Real(kind=8) :: b(m) = [ 2.5, 3.5, 1.5 ]
!  Real(kind=8) :: x(m) = 0.d0 
!Contains
!Subroutine Elimination ( )  !// ��˹��ȥ
!  Implicit none 
!  Integer :: i, j, k
!  Real(kind=8) :: mult, arr(m), arr_b
!
!  Write ( *,'(1x,a)' ) '������ȥǰ��������Ҷ���Ϊ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m ), b(i)
!  End do 
!  
!  Do j = 1, m - 1
!    !//---------------����Ԫ--------------
!    arr(:) = a(j,:); arr_b = b(j)
!    k = maxloc( abs( a(j:m,j) ), dim = 1 )
!    k = k + j - 1
!    a(j,:) = a(k,:); b(j) = b(k)
!    a(k,:) = arr(:); b(k) = arr_b
!    !//----------------------------------
!    Do i = j + 1, m
!      mult = a(i,j) / a(j,j)
!      Do k = j, m 
!        a(i,k) = a(i,k) - mult * a(j,k)
!      End do 
!      b(i) = b(i) - mult * b(j)
!    End do 
!  End do 
!
!  Write ( *,'(1x,a)' ) '������ȥ����������Ҷ���Ϊ��'
!  Do i = 1, m 
!    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m ), b(i)
!  End do 
!  
!End subroutine Elimination 
!
!Subroutine BackSubstitution ( )
!  Implicit none 
!  Integer :: i, j 
!  
!  Do i = m, 1, -1
!    Do j = i + 1, m 
!      b(i) = b(i) - a(i,j) * x(j)
!    End do 
!    x(i) = b(i) / a(i,i)
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
!Program GaussianElimination
!  Use mod 
!  Implicit none 
!  call Elimination ( )
!  call BackSubstitution ( )
!End program GaussianElimination