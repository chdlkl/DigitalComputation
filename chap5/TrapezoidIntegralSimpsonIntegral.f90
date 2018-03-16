Program main !// 辛普森积分
  Implicit none 
  Real(kind=8), parameter :: a = 0.d0, b = 1.d0   !//积分区间为[a,b], 被积函数为f(x) = x*e^x
  Integer, parameter :: n = 100  !// 区间等份
  Real(kind=8), parameter :: h = ( b - a ) / ( n*1.d0 )
  Real(kind=8) :: s = 0.d0  !// 积分结果
  Real(kind=8) :: x(0:n), tmp, y0, y1, y2  !// x:节点坐标
  Integer :: i
  
  Do i = 0, n 
    x(i) = a + i * h
  End do 
  
  !// 计算辛普森积分
  Do i = 1, n
    y0 = x(i-1) * exp( x(i-1) )
    tmp = ( x(i-1) + x(i) ) / 2.d0
    y1 = tmp * exp( tmp )
    y2 = x(i) * exp( x(i) )
    s = s + h/2.d0 * ( y0 + 4.d0 * y1 + y2 ) / 3.d0  !// 辛普森积分公式
  End do
  
  Write ( *,'(1x,A,g0)' ) '辛普森积分结果为: ', s
  
  !// 计算梯形积分
  s = 0.d0
  Do i = 1, n
    y0 = x(i-1) * exp( x(i-1) )
    y1 = x(i) * exp( x(i) )
    s = s + h * ( y0 + y1 ) / 2.d0  !// 梯形积分公式
  End do
  
  Write ( *,'(1x,A,g0)' ) '梯形积分结果为: ', s
  
  Write ( *,'(1x,A,g0)' ) '精确积分结果为: ', 1.d0
  
End program main