!// y'' = y' + cosy
!// y(0) = 0
!// y(pi) = 1
Module NoLinBVP
  Implicit none
  Integer, parameter :: n = 99  !// 待求点的个数
  Real(kind=8), parameter :: ya = 0.d0, yb = 1.d0  !// 边值
  Real(kind=8), parameter :: t0 = 0.d0, t1 = acos(-1.d0)
  Real(kind=8), parameter :: h = (t1-t0) / (n+1)
Contains
  Subroutine solveLinBVP()
    Implicit none
    Integer :: i
    Integer :: nloop = 20  !// 迭代次数
    Real(kind=8) :: y(n), deltay(n)
    
    y = 0.d0  !// 赋初值
    Do i = 1, nloop
      call Inv( y, deltay, n )
      y = y - deltay  !// 牛顿迭代法求解非线性方程组
    End do
    
    !// output
    open( 101, file = 'output.dat' )
    Do i = 1, n
      write( 100,* ) t0 + dble(i)*h, y(i)
    End do
    close( 101 )
  End subroutine solveNoLinBVP
  
  Subroutine Inv( y, deltay, n )
    Implicit none
    Integer :: i, j, k
    Real(kind=8) :: mult
    Real(kind=8) :: y(n), b(n), a(n,n), deltay(n)
    
    !// 系数矩阵a
    a = 0.d0  
    Do i = 1, n
      a(i,i) = -2.d0 + h*h*sind( y(i) )
    End do
    Do i = 1, n-1
      a(i,i+1) = 1.d0 - h / 2.d0
      a(i+1,i) = 1.d0 + h / 2.d0
    End do
    
    !// 右端项
    b = 0.d0
    b(1) = ( 1.d0 + h/2.d0 ) * ya - 2.d0 * y(1) + ( 1 - h/2.d0 ) * y(2) - h*h*cosd( y(1) )
    b(n) = ( 1.d0 + h/2.d0 ) * y(n-1) - 2.d0 * y(n) + ( 1 - h/2.d0 ) * yb - h*h*cosd( y(n) )
    Do i = 2:n-1
      b(i) = ( 1.d0 + h/2.d0 ) * y(i-1) - 2.d0 * y(i) + ( 1 - h/2.d0 ) * y(i+1) - h*h*cosd( y(i) )
    End do
    
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
   deltay = 0.d0
   Do i = n, 1, -1
     Do j = i + 1, n 
       b(i) = b(i) - a(i,j) * deltay(j)
     End do 
     deltay(i) = b(i) / a(i,i)
   End do
  End subroutine Inv
End module NoLinBVP
  
Program main
  use NoLinBVP
  call solveNoLinBVP
End program main