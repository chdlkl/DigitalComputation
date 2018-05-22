Program QR
  Implicit none
  Integer :: i, j
  Integer, parameter :: m = 4, n = 3
  Real(kind=8), parameter :: A(m,n) = reshape( [ 1.,2.,2.,-4.,3.,2.,2.,5.,-2.,6.,-4.,3. ],[m,n] )
  Real(kind=8) :: Q(m,n) = 0.d0, R(n,n) = 0.d0, y(m), A0(m,n)

  Do j = 1, n
    y = A(:,j)
    Do i = 1, j - 1
      R(i,j) = dot_product(Q(:,i),y)
      y = y - R(i,j)*Q(:,i)
    End do
    R(j,j) = sqrt( dot_product(y,y) )
    Q(:,j) = y / R(j,j)
  End do
  
  Write(*,'(1x,A)') 'The matrix Q is:'
  
  Do i = 1, m
    Write(*,'(*(f14.8))') Q(i,:)
  End do
  
  Write(*,'(1x,A)') 'The matrix R is:'

  Do i = 1, n
    Write(*,'(*(f14.8))') R(i,:)
  End do
  
  A0 = matmul(Q,R)
  Write(*,'(1x,A)') 'The matrix A is:'
  Do i = 1, m
    Write(*,'(*(f14.8))') A0(i,:)
  End do
  
End program QR