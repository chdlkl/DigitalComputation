Program QR
  Implicit none
  Integer :: i, j
  Integer, parameter :: m = 3, n = 2
  Real(kind=8), parameter :: A(m,n) = reshape( [ 1.,2.,2.,-4.,3.,2. ],[m,n] )
  Real(kind=8) :: Q(m,n) = 0.d0, R(n,n) = 0.d0
  Real(kind=8) :: y(m,1), qq(m,1), tmp(m,1), R1(1,1)
  Real(kind=8), external :: norm_2

  Do j = 1, n
    y(:,1) = A(:,j)
    Do i = 1, j - 1
      tmp(:,1) = A(:,j)
      R1 = matmul( transpose(qq),tmp )
      R(i,j) = R1(1,1)
      y = y - R(i,j)*qq
    End do
    R(j,j) = norm_2(y,m)
    qq = y / R(j,j)
    Q(:,j) = qq(:,1)
  End do
  
  Write(*,'(1x,A)') 'The matrix Q is:'
  Do i = 1, m
    Write(*,'(*(f12.5))') Q(i,:)
  End do
  Write(*,'(1x,A)') 'The matrix R is:'
  Do i = 1, n
    Write(*,'(*(f12.5))') R(i,:)
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
