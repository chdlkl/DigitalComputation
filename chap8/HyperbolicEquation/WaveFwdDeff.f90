!// Utt = c^2*Uxx   a<=x<=b, t>=0  c = 2
!// U(x,0) = f(x) = sin(pi*x)  a<=x<=b
!// Ut(x,0) = g(x) = 0  a<=x<=b
!// U(a,t) = L(t)  t>=0
!// U(b,t) = R(t)  t>=0
Module FwdDeff
  Implicit none
  Real(kind=8), parameter :: xL = 0.d0, xR = 1.d0, c = 2.d0
  Real(kind=8), parameter :: y0 = 0.d0, yT = 1.d0
  Real(kind=8), parameter :: pi = acos(-1.d0)
  Integer, parameter :: numH = 100, numT = 2000 !// �ռ䲽����ʱ�䲽��
Contains
  Subroutine solve()
    Implicit none
    Integer :: i, j
    Real(kind=8) :: h, k, x, t, sigma  !// DΪ��ɢϵ����hΪ�ռ䲽����kΪʱ�䲽��.  h��kҪ����: c*k <= h��Ϊ�����̵�CFL����
    Real(kind=8) :: u(0:numH,0:numT), Lside(0:numT), Rside(0:numT)
    
    h = ( xR-xL ) / numH
    k = ( yT-y0 ) / numT
    sigma = c*k / h
    If ( c*k >= h ) then
      write( *,'(1x,a)' ) ' ������CFL����, ������ֹ... '
      stop
    End if
    
    !// �߽�����
    Do i = 0, numT
      t = y0 + dble(i) * k
      Lside(i) = exp(-2.d0*t)
      Rside(i) = exp(-1.d0-2.d0*t)
    End do
    u(0,:) = Lside(0:numT)
    u(numH,:) = Rside(0:numT)
    
    !// ��ʼλ��f(x,0) = exp(-x)
    u(:,0) = exp( -1.d0*(xL + [0:numH]*h) )
    
    !// ���ó�ʼλ�Ƽ�����һ��ʱ�̵�λ��
    Do i = 1, numH-1
      x = xL + dble(i) * h
      u(i,1) = ( 1.d0 - sigma**2 ) * u(i,0) + k*g( x ) + sigma*sigma/2.d0 * ( u(i-1,0) + u(i+1,0) )
    End do
    
    Do j = 2, numT
      Do i = 1, numH-1
        u(i,j) = 2*( 1.d0 - sigma**2 ) * u(i,j-1) + sigma*sigma*( u(i-1,j-1) + u(i+1,j-1) ) - u(i,j-2)
      End do
    End do
    
    !// output
    open( 100, file = 'output.dat' )
    Do i = 0, numT
      t = y0 + dble(i) * k
      Do j = 0, numH
        x = xL + dble(j) * h
        write( 100,'(5f15.8)' ) t, x, u(j,i), exp( -x - 2.d0*t ), abs( exp( -x - 2.d0*t ) - u(j,i) ) / exp( -x - 2.d0*t ) * 100
      End do
    End do
    close( 100 )
    
  End subroutine solve
    
  Real(kind=8) function g( x )  !// ��ʼ�ٶ�
    Implicit none
    Real(kind=8), intent(in) :: x
    g = -2.d0*exp( -x )
  End function g
  
End module FwdDeff

Program main
  use FwdDeff
  Implicit none
  call solve
End program main  