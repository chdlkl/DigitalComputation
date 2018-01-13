Program Dichotomy  !// 二分法求解方程f(t) = cos(t) - t在区间[0,1]上的根
  Implicit none
  Integer :: i
  Integer, parameter :: n = 40 !// 迭代次数
  Real(kind=8) :: a = 0.d0, b = 1.d0, c, t
  Real(kind=8) :: fa, fb, fc

  fa = cos(a) - a   !// 方程为: f(x) = cos(x) - x
  fb = cos(b) - b
  If ( fa * fb > 0.d0 ) print*, 'error: f(a)*f(b) < 0 not satisfied!'

  Do i = 1, n 
    c = ( b + a ) / 2.d0
    fc = cos(c) - c
    If ( fc == 0.d0 ) exit
    If ( fa * fc > 0.d0 ) Then
      a = c; fa = fc
    Else If ( fa * fc < 0.d0 ) Then 
      b = c ; fb = fc
    End if 
  End do  
  t = ( b + a ) / 2.d0
  Write (*,'(1x,A,g0)') 't = ', t
  Write (*,'(" cos(t) - t = ",g0)') cos(t) - t
  
End program Dichotomy 