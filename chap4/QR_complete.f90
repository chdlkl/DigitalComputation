Subroutine calQR( R, m, n, Q ) 
  Implicit none
  Real(kind=8), external :: norm_2, sign0
  Integer :: m, n, i, j
  Real(kind=8) :: R(m,n), Q(m,m), I1(m,m), Q1(m,m), U(m,1), UTU(1,1)

  Q = 0.d0
  I1 = 0.d0
  Do i = 1, m
    Q(i,i) = 1.d0
    I1(i,i) = 1.d0
  End do

  Do j = 1, n
    U = 0.d0
    Do i = j, m
      U(i,1) = R(i,j)
    End do
    U(j,1) = U(j,1) + norm_2(U,m) * sign0( U(j,1) )
    UTU = matmul( transpose(U),U )
    Q1 = I1 - 2.d0 * matmul( U,transpose(U) ) / UTU(1,1)
    R = matmul( Q1,R )
    Q = matmul( Q1,Q )
    Q = transpose( Q )
  End do

End subroutine calQR

Real(kind=8) function norm_2( U,m )
  Implicit none
  Integer m, i
  Real(kind=8) :: U(m,1)
  
  norm_2 = 0.d0
  Do i = 1, m
    norm_2 = norm_2 + U(i,1)**2
  End do
  norm_2 = sqrt(norm_2)
  
End function norm_2

Real(kind=8) function sign0(t)
  Implicit none
  Real(kind=8) :: t
  Real(kind=8), parameter :: eps = 1.d-12
  
  If ( abs(t) < eps ) then
    sign0 = 1.d0
  else
    sign0 = t / abs(t)
  End if

End function sign0