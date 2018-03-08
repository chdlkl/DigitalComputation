Program QR
  Implicit none
  Integer :: i, j
  Integer, parameter :: m = 4, n = 3
  Real(kind=8), parameter :: A(m,n) = reshape( [ 1.,2.,2.,-4.,3.,2.,2.,5.,-2.,6.,-4.,3. ],[m,n] )
  Real(kind=8) :: Q(m,n) = 0.d0, R(n,n) = 0.d0
  Real(kind=8) :: y(m), qq(m), tmp(m)
  Real(kind=8), external :: norm_2

  Do j = 1, n
    y(:) = A(:,j)
    Do i = 1, j - 1
      tmp = A(:,j)
      qq = Q(:,i)
      R(i,j) = dot_product(qq,tmp)
      y = y - R(i,j)*qq
    End do
    R(j,j) = norm_2(y,m)
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
  
End program QR
  
Real(kind=8) function norm_2( y_tmp,m )
  Implicit none
  Integer :: i
  Integer, intent(in) :: m
  Real(kind=8), intent(in) :: y_tmp(m,1)

  norm_2 = 0.d0
  Do i = 1, m
    norm_2 = norm_2 + y_tmp(i,1)**2
  End do
  norm_2 = sqrt(norm_2)

End function norm_2