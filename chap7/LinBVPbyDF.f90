!// y'' = 4y
!// y(0) = 1
!// y(1) = 3
Module LinBVP
  Implicit none
  Integer, parameter :: n = 99  !// 待求点的个数
  Real(kind=8), parameter :: ya = 1.d0, yb = 3.d0  !// 边值
  Real(kind=8), parameter :: t0 = 0.d0, t1 = 1.d0
  Real(kind=8), parameter :: h = (t1-t0) / (n+1)
Contains
  Subroutine solveLinBVP()
    Implicit none
    Integer :: i, j, k
    Real(kind=8) :: y(n), b(n), a(n,n)
    Real(kind=8) :: mult
    
    !// 系数矩阵a
    a = 0.d0  
    Do i = 1, n
      a(i,i) = -4.d0*h*h - 2.d0
    End do
    Do i = 1, n-1
      a(i,i+1) = 1.d0
      a(i+1,i) = 1.d0
    End do
    
    !// 右端项
    b = 0.d0
    b(1) = -1.d0
    b(n) = -3.d0
    
    !// 解方程
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
     write( 100,* ) mult, y(i), ( 3.d0-exp(-2.d0) ) * exp(2.d0*mult) / ( exp(2.d0)-exp(-2.d0) ) + ( exp(2.d0)-3.d0 ) * exp(-2.d0*mult) / ( exp(2.d0)-exp(-2.d0) )
   End do
   close( 100 )
  End subroutine solveLinBVP
End module LinBVP
  
Program main
  use LinBVP
  call solveLinBVP
End program main