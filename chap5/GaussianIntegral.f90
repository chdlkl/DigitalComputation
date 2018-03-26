!// ע: �ο���Դhttp://fcode.cn/algorithm-73-1.html  
Module Gauss_Legendre !//��˹�����õ»��ָ�˹�㼰Ȩ�ص����ģ��
  Implicit none
  Integer, parameter :: n  = 11                          !// ��������˹��ĸ���
  Integer, parameter :: DP = selected_real_kind( p=13 )  !// ����kind����ֵ
  Real(kind=DP), parameter :: eps = 1.0e-15_DP           !// ��������
Contains
  Real(Kind=DP) function N_Legendre(x) !// ����n�����õ¶���ʽ
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a(n), x
    a(1) = x !// 1�����õ¶���ʽ
    a(2) = 1.5_DP*x*x - 0.5_DP !// 2�����õ¶���ʽ
    Do i = 3, n
      a(i) = ( dble(i+i-1)*x*a(i-1) - dble(i-1)*a(i-2) ) / dble(i) !// ���õ��ƹ�ϵ����n�����õ¶���ʽ
    End do
    N_Legendre=a(n) !//���ɵ�n�����õ¶���ʽ
  End function N_Legendre

  Real(Kind=DP) Function N1_Legendre(x)  !// ����n-1�����õ¶���ʽ 
    Implicit none
    Integer :: i
    Real (Kind=DP) :: a(n), x
    a(1) = x
    a(2) = 1.5_DP*x**2 - 0.5_DP
    Do i = 3, n - 1
      a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
    End Do
    N1_Legendre = a(n-1)     
  End function N1_Legendre
  
  Real(Kind=DP) function DN_Legendre(x)  !// ����n�����õ¶���ʽ�ĵ������ʽ
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a(n), x
    a(1) = x  !// 1�����õ¶���ʽ
    a(2) = 1.5_DP*x*x - 0.5_DP !// 2�����õ¶���ʽ
    Do i = 3, n
      a(i) = ( dble(i+i-1)*x*a(i-1) - dble(i-1)*a(i-2) ) / dble(i) !// ���õ��ƹ�ϵ����n�����õ¶���ʽ
    End Do
    DN_Legendre = ( a(n-1) - x*a(n) )*dble(n) / (1.0_DP - x*x ) 
  End function DN_Legendre

  Real(Kind=DP) function NewtonIteration(a, b) !// ţ�ٷ���⺯���Ľ�
    Implicit none
    Integer :: i
    Real(Kind=DP) :: a, b, x, xtmp
    Integer, parameter :: nloop = 2000
    !// a,b�Ǵ��ݽ����Ļ��ֺõ���һ������ڵ�����
    x = ( a + b ) / 2.d0  !// ��ʼ����ֵ
    i = 0
    Do 
      xtmp = x - N_Legendre(x) / DN_Legendre(x)   !// X(i+1) = Xi - f(Xi) / f'(Xi)  i = 1,2,...N
      i = i + 1
      If ( abs( xtmp-x ) < eps .and. i > nloop ) exit
      x = xtmp
    End do 
    NewtonIteration = x
  End function NewtonIteration

  Subroutine root_coeff ( f_root, f_coeff )  !// ����N�����õ¶���ʽ�ĸ���ȥ��Ȩ��ϵ��
    Implicit none
    Real(Kind=DP) :: m, nstep, f_root(n), f_coeff(n) !// ��������,��Сn��module��ʼ������
    Integer :: i, j
    Real(kind=DP), parameter :: h = 1.d-6
    j = 0   !// ��ֵ����ѭ�������ĳ�ֵ           
    m = -1.d0 + h   !// ���ü�����[-1��1] �����ޣ�������-1 
    nstep = nint(2.d0/h)
    Do i = 1, nstep   !// ���ѭ������Ӧ�����ɲ���0.000001�� ��,���㷽����2000000=2/0.000001     
      If ( N_Legendre(m)*N_Legendre(m+h) < 0 ) then   !// �����޴���ʼ�������ۼ�
        j = j + 1    !// ��¼���ǵڼ�����
        f_root(j) = NewtonIteration( m, m+h )!// ����ţ�ٷ��������ڷֺõ�һС������⣬����洢��fn��j��
        f_coeff(j) = 2.0_DP / ( dble(n) * N1_Legendre(f_root(j)) * DN_Legendre(f_root(j)) ) !// ���ù�ʽ�����˹���Ȩ��
        write (*,'(1x,a,g0)') '��˹�����: ', j
        write (*,'(1x,a,g0,2x,a,g0)') '��˹��: ', f_root(j), '��˹��Ȩ��', f_coeff(j)
        write(*,'(1x,a)') '------------------------------------------------------'
      End if
      m = m + h !// ִ����һ���ж�m��ǰ�ƽ�һ��
    End Do
  End subroutine root_coeff
  
  Real(Kind=DP) function func(x) !// ��������
    Implicit none
    Real(Kind=DP) :: x
    func = exp( -x*x/2.0_DP ) !// ÿ�μ����������޸ı�����������
  End function func

End module Gauss_Legendre

Program GaussianIntegral
  use Gauss_Legendre
  Implicit none
  Real (Kind=DP) :: f_root(n), f_coeff(n), x, a, b, answer
  Integer :: i
  Call root_coeff ( f_root, f_coeff ) !// �������˹����Ȩ�ص��Ӻ���
  a = -1.0_DP !// ��������
  b = 1.0_DP !// ��������   
  answer = 0.d0 !// ����ֽ������ʼֵ
  !// һ������[a,b]�ϵĻ��ֹ�ʽ
  !// Integral[f(x),a,b] = Integral[f( ((b-a)*t+b+a)/2 )] * ( b-a )/2. tΪN�����õ¶���ʽ�ĸ�
  Do i = 1, n
    answer = answer + f_coeff(i) * func( (a+b) / 2.0_DP + (b-a) / 2.0_DP * f_root(i) ) !// ��˹���õ�����ֹ�ʽ     
  End Do
  answer = answer * (b-a) / 2.0_DP
  !// ��ȷ��Ϊ1.71124878378430
  !// ��ֵ��Ϊ1.711248783784299  
  Write(*,'(1x,a,g0)') '��˹-���õ�����ֽ��: ', answer

End program GaussianIntegral
