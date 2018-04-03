!// y' = y + 8*y^2 - 9*y^3
!// y(0) = 0.5
!// t = [0,3]
Module BWDeuler
  Implicit none
Contains
  Subroutine solve()
    Implicit none
    Real(kind=8), parameter :: t0 = 0.d0, t1 = 3.d0
    Real(kind=8), parameter :: h = 0.01d0
    Real(kind=8), parameter :: n = nint( (t1-t0)/h ), eps = 1d-6
    Integer, parameter :: nloop = 1000
    Integer :: i, j
    Real(kind=8) :: y, y_tmp, t, y0
    
    y0 = 0.5d0
    open( 100, file = 'BWDeuler.dat' )
    write( 100,* ) t0, y0
    Do i = 1, n
      j = 0
      y_tmp = y0  !// 每一个节点初值的选取都是上一个节点的值，这样可以加快收敛速度
      Do
        y = y_tmp - func( h, y_tmp, y0 ) / DerivativeFunc( h, y_tmp )  !// 使用牛顿迭代法求解每一个节点的值
        j = j + 1
        If ( abs(y-y_tmp) < eps .or. j > nloop ) exit
        y_tmp = y
      End do
      y0 = y
      t = t0 + i*h
      write( 100,* ) t, y0
    End do
    close( 100 )
  End subroutine solve
  
  Real(kind=8) function func( h, y_tmp, y0 )  
    Implicit none 
    Real(kind=8), intent(in) :: h, y_tmp, y0 
    func = 9.d0 * h * y_tmp**3 - 8.d0 * h * y_tmp**2 + ( 1.d0 - h ) * y_tmp - y0
  End function func
  
  Real(kind=8) function DerivativeFunc( h, y_tmp )  
    Implicit none 
    Real(kind=8), intent(in) :: h, y_tmp
    DerivativeFunc = 27.d0 * h * y_tmp**2 - 16.d0 * h * y_tmp +  1.d0 - h 
  End function DerivativeFunc
End module BWDeuler 
  
  
Program main
  use BWDeuler
  call solve
End program main