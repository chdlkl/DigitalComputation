Program main !// �������λ��ֺ�����ɭ����
  Implicit none 
  Real(kind=8), parameter :: a = 0.d0, b = 1.d0   !//��������Ϊ[a,b], ��������Ϊf(x) = x*e^x
  Integer, parameter :: n = 100  !// ����ȷ�
  Real(kind=8), parameter :: h = ( b - a ) / ( n*1.d0 )
  Real(kind=8) :: s = 0.d0  !// ���ֽ��
  Real(kind=8) :: x(0:n), y(n-1), y1(n), y2(n-1), tmp, y0, ym, y2m  !// x:�ڵ�����
  Integer :: i, j
  
  Do i = 0, n 
    x(i) = a + i * h
    If ( i > 0  .and. i < n ) then
      y(i) = x(i) * exp( x(i) )
    End if
  End do 
  
  s = 0.d0
  !// ���㸴�����λ���
  y0 = x(0) * exp( x(0) )
  ym = x(n) * exp( x(n) )
  s = h * ( y0 + ym + 2.d0*sum(y) ) / 2.d0  !// �������λ��ֹ�ʽ
  Write ( *,'(1x,A,g0)' ) '�������λ��ֽ��Ϊ: ', s
  
  !// ���㸴������ɭ����
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
  !// ��ʱ��������ɭ���ֹ�ʽ�Ĳ����Ǹ������λ��ֵ�һ�룬����Ϊh/2.d0
  s = h / 2.d0 * ( y0 + y2m + 4.d0*sum(y1) + 2.d0*sum(y2) ) / 3.d0  !// ��������ɭ���ֹ�ʽ
  
  Write ( *,'(1x,A,g0)' ) '����ɭ���ֽ��Ϊ: ', s
  
  Write ( *,'(1x,A,g0)' ) '��ȷ���ֽ��Ϊ: ', 1.d0
  
End program main