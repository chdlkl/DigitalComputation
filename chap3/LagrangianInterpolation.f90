Program LagrangianInterpolation
  Implicit none
  Integer :: i, fileid 
  Integer, parameter :: n1 = 51
  Real(kind=8), parameter :: t0 = 0.d0, t1 = 1.d0
  Real(kind=8) :: x(n1), y(n1)
  
  !// 获取已知点数据
  Open ( newunit = fileid, file = 'origin.dat' )
  Do i = 1, n1
    x(i) = ( i-1 ) * ( t1-t0 ) / ( n1-1 )
    y(i) = exp( x(i) )
    Write ( fileid,'(2es)' ) x(i), y(i)
  End do 
  Close ( fileid )
  !// 求插值
  call GetInterpolation( x, y, n1, t0, t1 )
End program LagrangianInterpolation
  
Subroutine GetInterpolation( x, y, n1, t0, t1 )
  Implicit none
  Integer :: i, j, k, fileid
  Integer, intent(in) :: n1
  Integer, parameter :: n = 101
  Real(kind=8), intent(in) :: x(n1), y(n1), t0, t1 
  Real(kind=8) :: xx(n), yy(n), L1(n1), L2(n1), tmp
  
  Open ( newunit = fileid, file = 'chazhi.dat' )
  Do k = 1, n
    xx(k) = ( k-1 ) * ( t1-t0 ) / ( n-1 )
    Do i = 1, n1
      L1(i) = 1.d0
      Do j = 1, n1
        If ( j /= i ) then
          L1(i) = ( x(i) - x(j) ) * L1(i)  !// 求分母
        End if 
      End do
      
      L2(i) = 1.d0
      Do j = 1, n1
        If ( j /= i ) then
          L2(i) = ( xx(k) - x(j) ) * L2(i)  !// 求分子
        End if
      End do
    End do
    
    tmp = 0.d0
    Do i = 1, n1
      tmp = y(i) * L2(i) / L1(i) + tmp  !// 求插值项
    End do
    yy(k) = tmp
    Write ( fileid,'(2es)' ) xx(k), yy(k)
  End do
  Close ( fileid )
End subroutine GetInterpolation