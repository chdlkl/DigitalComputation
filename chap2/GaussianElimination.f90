!//---------------原方程-------------
!//     x + y + z + w = 10
!//    2x + 3y + z + w = 15
!//    3x - y + 2z - w = 3
!//    4x + y -3z + 2w = 5
!//----------------------------------
Module mod 
  Implicit none 
  Integer, parameter :: m = 4
  Real(kind=8) :: a(m,m) = [ 1.d0, 2.d0, 3.d0, 4.d0, 1.d0, 3.d0, -1.d0, 1.d0, 1.d0, 1.d0, 2.d0, -3.d0, 1.d0, 1.d0, -1.d0, 2.d0 ]
  Real(kind=8) :: b(m) = [ 10.d0, 15.d0, 3.d0, 5.d0 ]
  Real(kind=8) :: x(m) = 0.d0 
Contains
Subroutine Elimination ( )  !// 高斯消去
  Implicit none 
  Integer :: i, j, k 
  Real(kind=8), parameter :: eps = 1.d-8
  Real(kind=8) :: mult

  Write ( *,'(1x,a)' ) '经过消去前左端项与右端项为：'
  Do i = 1, m 
    Write ( *,'(*(f12.5))' ) ( a(i,j), j = 1, m ), b(i)
  End do 
  
  Do j = 1, m - 1
    If ( abs(a(j,j)) < eps ) Then
      Write ( *,'(1x,a)' ) ' The pivot is zero!'
      stop 
    End if 
    Do i = j + 1, m
      mult = a(i,j) / a(j,j)
      Do k = j, m 
        a(i,k) = a(i,k) - mult * a(j,k)
      End do 
      b(i) = b(i) - mult * b(j)
    End do 
  End do 

  Write ( *,'(1x,a)' ) '经过消去后左端项与右端项为：'
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

  Write ( *,'(1x,a)' ) '原方程解为：'
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