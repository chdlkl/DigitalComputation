Program euler
  Implicit none
  call euler_slu()
End program euler

Subroutine euler_slu()
  Implicit none
  Real(kind=8), external :: derfunc
  Integer, parameter :: n = 10 !// 分成100等份
  Real(kind=8) :: a = 0.d0, b = 1.d0  !// 时间区间
  Real(kind=8) :: y0, t, y, h
  Integer :: i

  y0 = 1.d0  !// 初值
  h = ( b - a ) / n
  open ( 100, file = 'euler_slu.dat' )
  write( 100,* ) a, y0
  Do i = 1, n
    t = a + h * ( i - 1 )
    y = y0 + h * derfunc( t, y0 ) 
    write( 100,* ) t+h, y
    y0 = y
  End do
  close( 100 )

  open ( 100, file = 'slu.dat' )
  Do i = 1, n+1
    t = a + h * ( i - 1 )
    y = 3.d0 * exp(t**2/2.d0) - t**2 - 2.d0
    write( 100,* ) t, y
  End do
  close( 100 )
        
End subroutine euler_slu
  
Real(kind=8) function derfunc( t, y )
  Implicit none
  Real(kind=8), intent(in) :: t, y
  derfunc = t * y + t**3
End function derfunc

!// 初值问题如下
!// y' = t*y + t^3
!// y(0) = y0
!// t = [0,1]
!// 解析解为：y(t) = 3*exp(t^2/2) - t^2 - 2
