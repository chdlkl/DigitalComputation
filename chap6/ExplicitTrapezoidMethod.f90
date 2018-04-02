Program ExplicitTrapezoidMethod
  Implicit none
  call TrapezoidMethod_slu()
End program ExplicitTrapezoidMethod

Subroutine TrapezoidMethod_slu()
  Implicit none
  Real(kind=8), external :: derfunc
  Integer, parameter :: n = 100 !// �ֳ�100�ȷ�
  Real(kind=8) :: a = 0.d0, b = 1.d0  !// ʱ������
  Real(kind=8) :: y0, t, y, h, y_tmp
  Integer :: i

  y0 = 1.d0  !// ��ֵ
  h = ( b - a ) / n
  open ( 100, file = 'TrapezoidMethod_slu.dat' )
  write( 100,* ) a, y0
  Do i = 1, n
    t = a + h * ( i - 1 )
    !// Wi+1 = Wi + h * ( f(Ti,Wi) + f( Ti+h,Wi+h*f(Ti,Wi) ) )
    y_tmp = derfunc( t, y0 )
    y = y0 + h * ( derfunc( t, y0 ) + derfunc( t + h, y0 + h*derfunc( t, y0 ) ) ) / 2.d0
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
        
End subroutine TrapezoidMethod_slu
  
Real(kind=8) function derfunc( t, y )
  Implicit none
  Real(kind=8), intent(in) :: t, y
  derfunc = t * y + t**3
End function derfunc

!// ��ֵ��������
!// y' = t*y + t^3; f = t*y + t^3
!// y(0) = y0
!// t = [0,1]
!// ������Ϊ��y(t) = 3*exp(t^2/2) - t^2 - 2
