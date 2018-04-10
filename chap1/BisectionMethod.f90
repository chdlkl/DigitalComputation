Program BisectionMethod  !// 二分法求解方程f(t) = cos(t) - t在区间[0,1]上的根
  Implicit none
  Real(kind=8), parameter :: eps = 1d-12 !// 精度控制
  Real(kind=8) :: a = 0.d0, b = 1.d0, c, t
  Real(kind=8) :: fa, fb, fc

  fa = cos(a) - a   !// 方程为: f(x) = cos(x) - x
  fb = cos(b) - b
  If ( fa * fb > 0.d0 ) then 
    print*, 'error: f(a)*f(b) < 0 not satisfied!'
    stop
  end if

  Do while ( (b - a) > eps )  
    c = ( b + a ) / 2.d0
    fc = cos(c) - c
    If ( fc == 0.d0 ) exit
    If ( fa * fc > 0.d0 ) then
      a = c; fa = fc
    else
      b = c ; fb = fc
    End if 
  End do  
  t = ( b + a ) / 2.d0
  Write (*,'(1x,A,f20.15)') 't =', t
  Write (*,'(" cos(t) - t =",f20.15)') cos(t) - t
  
End program BisectionMethod 