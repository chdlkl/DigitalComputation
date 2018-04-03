!---------------------------------------------------------
!             �˳����õ���AdamsBashforth(�ಽ��巨)���
!                       y'=-3*y
!                       y(0) = 1
!                       t = [0,2]
!                   �����⣺y = e^(-3t)
!---------------------------------------------------------
Module AdamsMoulton
  Implicit none
Contains
  Subroutine solve()  
    Real(kind=8), parameter :: T0 = 0.d0, T1 = 2.d0 !//nΪ������T1��T0Ϊ��������
    Real(kind=8), parameter :: h = 1d-2  !//hΪ����
    Integer, parameter :: n = 3
    Integer, parameter :: num = nint((T1-T0)/h)  !//numΪ�ܵ���
    Real(kind=8) :: a(0:n), b(0:n), sum, t, sum1, b1, jL, u1(0:num), t11, t22, u11(0:num), tmp
    !//����a������ɺ������ĵ���ϵ��������b������ձ��ʽǰ��ϵ��
    !//sum,sum1,b1��Ϊ��ʼ�������������ѭ��ֵ��t11,t22Ϊ��ʱ���������е�������
    !//u1���������ֵ�⣬u11�������巨�������ֵ��
    Integer :: i, j, L
  
    write(*,*) "----------���Adams�ڲ幫ʽϵ��-------------"
    a(0) = 1.d0
    Do i = 1, n
      sum = 0.d0; t = 2.d0
      Do j = i, 1, -1
        sum = sum + a(j-1) / (t)
        t = t + 1.d0
      End Do
      a(i) = 0.d0 - sum  !//�˴�����巨��һ����a(i)=0.-sum;����巨Ϊa(i)=1.-sum
    End Do
    write(*,*) a
  
    Print*,"----------��֤Adams�ڲ幫ʽϵ��-------------"
    Print*
    sum1 = 0.d0  !//��֤Adams�ڲ幫ʽϵ��ֵ�ֽ��Ƿ���ȷ
    Do i = 0, n
      sum1 = sum1 + a(n-i) / dble(1+i)
    End do
    If( abs(sum1-0.0) < 1e-6 ) then   !//sum1=0.0,˵��Adams�ڲ幫ʽϵ���ֽ���ȷ
      write(*,*) "----------Adams�ڲ幫ʽϵ���ֽ���ȷ---------"  
    else
      write(*,*) "---------Adams��幫ʽϵ���ֽⲻ��ȷ--------"
      stop
    End if
  
    Print*, '--------------------------------------------'
    Write(*,*) "---------------���ϵ��Bklֵ----------------"
    Print*, '--------------------------------------------'
  
    Do L = 0, n
      b1 = 0.d0
      Do j = L, n
        call calculate_jL( j, L, jL )  !//���(s(s-1)...(s-j+1)/j!)
        b1 = b1 + (-1)**L * a(j) * jL
      End do
      b(L) = b1
    End do
    write(*,*) b*24  !//24��Ϊ�˽�ϵ����ʾ��������û�������ر���  

    open( 100, file = 'waicha.dat', status = 'old')  !//��ȡAdams��巨�������ֵ
    Do i = 0, num
      read(100,*) tmp, u11(i)
      If( i<n ) u1(i) = u11(i)   !//ǰ3����ֵ������㷨����
    End do
    close( 100 )
  
    !//AdamsMoulton��ʽ
    Do i = n, num  !���ڱ������õ���3�ף������ڲ巨��ǰ3����ֵҪ����巨�������ڲ巨Ҫ����ǰ4����ֵ��
      t11 = u11(i)
      Do
        t22 = u1(i-1) + h * ( b(0)*(-3.*t11) + b(1)*(-3.*u1(i-1)) + b(2)*(-3.*u1(i-2)) + b(3)*(-3.*u1(i-3)) )
        If( abs(t22-t11) < 1d-6 ) then  
          exit  !//������������Ļ�������ѭ������ֵ��u1(i)
        else
          t11 = t22  !//���������Ҫ�󣬼�������
        End if
      End do
      u1(i) = t22
    End do
  
    open( 100, file = 'neicha.dat', status='unknown')
    Do i = 0, num
      write( 100,'(4f12.8)' ) i*h, u1(i), exp(-3.d0*h*dble(i)), abs( (u1(i)-exp(-3.d0*h*dble(i))) / exp(-3.d0*h*dble(i)) ) * 100
    End do
    close(100)
  End subroutine solve
  
  Subroutine calculate_jL( j, L, jL )
    Implicit None
    Integer :: i
    Integer, intent(in) :: j, L
    Real(kind=8), intent(inout) :: jL
    Real(kind=8) :: Product1, Product2
 
    If( L == 0 ) then
      jL = 1.d0
      return
    else
      Product1 = 1.d0; Product2 = 1.d0
      Do i = 1, L
        Product1 = Product1 * i  !//������s(s-1)...(s-j+1)
      End do
      Do i = j, j-L+1, -1
        Product2 = Product2 * i  !//���j�Ľ׳�
      End do
      jL = Product2 / Product1
    End If
  End subroutine calculate_jL
End module AdamsMoulton

Program main
  use AdamsMoulton
  Implicit none
  call solve()
End program main
  
