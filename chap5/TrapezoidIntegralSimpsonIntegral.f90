Program main !// ����ɭ����
  Implicit none 
  Real(kind=8), parameter :: a = 0.d0, b = 1.d0   !//��������Ϊ[a,b], ��������Ϊf(x) = x*e^x
  Integer, parameter :: n = 100  !// ����ȷ�
  Real(kind=8), parameter :: h = ( b - a ) / ( n*1.d0 )
  Real(kind=8) :: s = 0.d0  !// ���ֽ��
  Real(kind=8) :: x(0:n), tmp, y0, y1, y2  !// x:�ڵ�����
  Integer :: i
  
  Do i = 0, n 
    x(i) = a + i * h
  End do 
  
  !// ��������ɭ����
  Do i = 1, n
    y0 = x(i-1) * exp( x(i-1) )
    tmp = ( x(i-1) + x(i) ) / 2.d0
    y1 = tmp * exp( tmp )
    y2 = x(i) * exp( x(i) )
    s = s + h/2.d0 * ( y0 + 4.d0 * y1 + y2 ) / 3.d0  !// ����ɭ���ֹ�ʽ
  End do
  
  Write ( *,'(1x,A,g0)' ) '����ɭ���ֽ��Ϊ: ', s
  
  !// �������λ���
  s = 0.d0
  Do i = 1, n
    y0 = x(i-1) * exp( x(i-1) )
    y1 = x(i) * exp( x(i) )
    s = s + h * ( y0 + y1 ) / 2.d0  !// ���λ��ֹ�ʽ
  End do
  
  Write ( *,'(1x,A,g0)' ) '���λ��ֽ��Ϊ: ', s
  
  Write ( *,'(1x,A,g0)' ) '��ȷ���ֽ��Ϊ: ', 1.d0
  
End program main