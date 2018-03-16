Program main !// 复合梯形积分和辛普森积分
  Implicit none 
  Real(kind=8), parameter :: a = 0.d0, b = 1.d0   !//积分区间为[a,b], 被积函数为f(x) = x*e^x
  Integer, parameter :: n = 100  !// 区间等份
  Real(kind=8), parameter :: h = ( b - a ) / ( n*1.d0 )
  Real(kind=8) :: s = 0.d0  !// 积分结果
  Real(kind=8) :: x(0:n), y(n-1), y1(n), y2(n-1), tmp, y0, ym, y2m  !// x:节点坐标
  Integer :: i, j
  
  Do i = 0, n 
    x(i) = a + i * h
    If ( i > 0  .and. i < n ) then
      y(i) = x(i) * exp( x(i) )
    End if
  End do 
  
  s = 0.d0
  !// 计算复合梯形积分
  y0 = x(0) * exp( x(0) )
  ym = x(n) * exp( x(n) )
  s = h * ( y0 + ym + 2.d0*sum(y) ) / 2.d0  !// 复合梯形积分公式
  Write ( *,'(1x,A,g0)' ) '复合梯形积分结果为: ', s
  
  !// 计算复合辛普森积分
  s = 0.d0
  i = 0
  Do j = 1, 2*n
    If ( mod(j,2) == 1 ) then
      i = i + 1
      tmp = dble(j)*h / 2.d0
      y1(i) = tmp * exp( tmp )
    End if
  End do
    
  i = 0
  Do j = 1, 2*n - 1
    If ( mod(j,2) == 0 ) then
      i = i + 1
      tmp = dble(j)*h / 2.d0
      y2(i) = tmp * exp( tmp )
    End if
  End do

  y0 = x(0) * exp( x(0) )
  y2m = x(n) * exp( x(n) )
  !// 此时复合辛普森积分公式的步长是复合梯形积分的一半，所以为h/2.d0
  s = h / 2.d0 * ( y0 + y2m + 4.d0*sum(y1) + 2.d0*sum(y2) ) / 3.d0  !// 复合辛普森积分公式
  
  Write ( *,'(1x,A,g0)' ) '辛普森积分结果为: ', s
  
  Write ( *,'(1x,A,g0)' ) '精确积分结果为: ', 1.d0
  
End program main