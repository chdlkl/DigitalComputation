Program ConjugateGradientMethods
  Implicit none
  Integer :: i
  Integer, parameter :: n = 3
  Real(kind=8) :: A(n,n) = [ 4.d0, -2.d0, 2.d0, -2.d0, 2.d0, -4.d0, 2.d0, -4.d0, 11.d0 ]
  Real(kind=8) :: b(n,1) = [ 6.d0,-10.d0,27.d0 ]
  Real(kind=8) :: alpha = 0.d0, beta = 0.d0
  Real(kind=8) :: x(n,1) = 0.d0, d(n,1) = 0.d0, r(n,1) = 0.d0, rtmp(n,1) = 0.d0, tmp(1,1) = 0.d0
  
  d = b - matmul( A,x )
  r = d

  Do i = 0, n - 1
    If ( all(r<0.0) ) then
      exit
    End if
    rtmp = r
    tmp = matmul( transpose(r),r ) / matmul( matmul( transpose(d),A ), d )
    alpha = tmp(1,1)
    x = x + alpha * d
    r = r - alpha * matmul( A,d )
    tmp = matmul( transpose(r),r ) / matmul( transpose(rtmp),rtmp )
    beta = tmp(1,1)
    d = r + beta * d
  End do
  
  Do i = 1, n
    Write ( *,'(f12.5)' ) x(i,1)
  End do
End program ConjugateGradientMethods