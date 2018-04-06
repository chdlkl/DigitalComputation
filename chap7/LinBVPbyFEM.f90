!// y'' = 4y
!// y(0) = 1
!// y(1) = e^2
Module LinBVPbyFEM
  Implicit none
  Real(kind=8), parameter :: t0 = 0.d0, t1= 1.d0
  Real(kind=8), parameter :: ya = 1.d0, yb = exp(2.d0)
  Integer, parameter :: n = 99
  Real(kind=8), parameter :: h = ( t1-t0 ) / ( n+1 )
Contains
  Subroutine solveByFEM()
    Implicit none
    integer :: i, j, k
    Real(kind=8) :: a(n,n), b(n), y(n)
    Real(kind=8) :: alpha, beta, mult
    
    !// 系数矩阵a
    alpha = 8.d0*h/3.d0 + 2.d0/h
    beta = 2.d0*h/3.d0 - 1.d0/h
    a = 0.d0
    Do i = 1, n
      a(i,i) = alpha
    End do
    Do i = 1, n-1
      a(i,i+1) = beta
      a(i+1,i) = beta
    End do
    
    !// 右端项
    b = 0.d0
    b(1) = -1.d0*ya*beta
    b(n) = -1.d0*yb*beta
    
    !// 用gauss消去法求解
    Do j = 1, n - 1
      Do i = j + 1, n
        mult = a(i,j) / a(j,j)
        Do k = j, n 
          a(i,k) = a(i,k) - mult * a(j,k)
        End do 
        b(i) = b(i) - mult * b(j)
      End do 
    End do 
    
   !// 回代
   y = 0.d0
   Do i = n, 1, -1
     Do j = i + 1, n 
       b(i) = b(i) - a(i,j) * y(j)
     End do 
     y(i) = b(i) / a(i,i)
   End do
   
   !// output
   open( 100, file = 'output.dat' )
   Do i = 1, n
     mult = t0 + dble(i)*h
     write( 100,* ) mult, y(i), exp(2.d0*mult)  !// 解析解为: y(t) = exp(2t)
   End do
   close( 100 )
  End subroutine solveByFEM
End module LinBVPbyFEM
  
Program main
  use LinBVPbyFEM
  call solveByFEM
End program main