Module mod
  Implicit none
  Integer, parameter, public :: m = 4, n = 2
  !Real(kind=8), public :: A(m,n) = reshape( [1.,2.,2.,-4,-4.,3.,2.,-2,1.,2.,3.,4.],[m,n] )
  Real(kind=8), public :: A(m,n) = reshape( [1.,2.,2.,-4.,-4.,3.,2.,-2.],[m,n] )
Contains
  Subroutine calHouseholderMatrix
    Implicit none
    Integer :: i, j, k, length
    Real(kind=8) :: norm_x
    Real(kind=8), allocatable :: w(:), v(:,:), x(:)
    Real(kind=8), allocatable :: P(:,:), tmpH(:,:), tmpI(:,:)
    Real(kind=8) :: Q(m,m), R(m,n), II(m,m)
    
    R = A
    k = 0
    Do i = 1, n
      length = m - k
      allocate( x(length), v(length,1), w(length) )
      allocate( P(length,length), tmpH(length,length), tmpI(length,length) )
      
      x = R(i:m,i)
      w = 0.d0
      
      norm_x = dot_product( x,x )
      w(1) = sqrt( norm_x )  !// w = [ |x|2, 0, ..., 0 ]
      v(:,1) = w - x  !// v = w - x
      P = matmul( v,transpose(v) ) / dot_product( v(:,1),v(:,1) )  !// P = v*v'/(v'*v)
      
      tmpI = 0.d0
      Do j = 1, length
        tmpI(j,j) = 1.d0
      End do
      tmpH = tmpI - 2.d0 * P  !// H = I - 2P
      
      If ( i < 2 ) then
        R = matmul( tmpH, A )
        Q = tmpH
      else
        II = 0.d0
        Do j = 1, m
          II(j,j) = 1.d0
        End do
        II(i:m,i:m) = tmpH
        Q = matmul( Q, II )
        R = matmul( II, R )
      End if
      
      k = k + 1
      Deallocate( x, v, w, P, tmpH, tmpI )
    End do
    
    Write ( *,'(1x,a)' ) 'The matrix Q is:'
    Do i = 1, m
      Write (*,'(*(f14.8))') Q(i,:)
    End do
    
    Write ( *,'(1x,a)' ) 'The matrix R is:'
    Do i = 1, m
      Write (*,'(*(f14.8))') R(i,:)
    End do
    
    Write ( *,'(1x,a)' ) 'The matrix A is:'
    A = matmul( Q,R )
    Do i = 1, m
      Write (*,'(*(f14.8))') A(i,:)  !// A = QR
    End do
    
    Q = matmul(Q,transpose(Q))
    Write ( *,'(1x,a)' ) 'The matrix Q*Q'' is:'
    Do i = 1, m
      Write (*,'(*(f14.8))') Q(i,:)  !// Q具有正交性。Q*Q' = I
    End do
    
  End subroutine calHouseholderMatrix
End module mod
  
Program HouseholderMatrix
  use mod
  Implicit none
  call calHouseholderMatrix
End program HouseholderMatrix