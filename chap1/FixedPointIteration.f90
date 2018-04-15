Program FixedPointIteration  !// ʹ�ò��������������ⷽ�� x^3 + x - 1 = 0 �����
  Implicit none 
  Real(kind=8), parameter :: eps = 1.d-8
  Real(kind=8), external :: func
  Real(kind=8) :: xtemp, x0 = 0.5d0  !// x0Ϊ������ʼֵ
  
  Do 
    xtemp = func ( x0 )
    If ( abs( xtemp - x0 ) < eps ) exit 
    x0 = xtemp
  End do 
  Write ( *,'(a,f13.9)' ) ' the result of fixed point iteration is', xtemp

End program FixedPointIteration
  
Real(kind=8) function func ( x )
  Implicit none 
  Real(kind=8), intent(in) :: x 
  func = ( 1.d0 + 2.d0 * x * x * x ) / ( 1.d0 + 3.d0 * x * x )
End function func 
!// ���� x^3 + x - 1 = 0 ��Ҫ�Ĳ������ʽ������
!// ��һ�֣�x = 1 - x^3 ��һ�ָ�ʽ������
!// �ڶ��֣�x = ( 1 - x^3 )^(1/3)�����ָ�ʽ����
!// �����֣����Ǵ˴����е�������ʽ��x = ( 1 + 2*x^3 ) / ( 1 + 3*x^2 )�����ָ�ʽ�������������������ٶȿ죬���ȸ�