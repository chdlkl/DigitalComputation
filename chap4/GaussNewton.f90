!// ������ԭ���Լ�������������wiki: https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm#Example
Module varGS
  Implicit None
  !//NparasΪ����������ndataΪ���ݸ�����NitersΪ����������
  Integer, parameter :: Nparas = 2, ndata = 7, Niters = 50, Fileid = 101
  !//xobsΪ�۲������Ա�����yobsΪ�۲���������������߾�Ϊ��֪����
  Real(kind=8), parameter :: xobs(ndata) = [0.038d0, 0.194d0, 0.425d0, 0.626d0, 1.253d0, 2.500d0, 3.740d0]
  Real(kind=8), parameter :: yobs(ndata) = [0.050d0, 0.127d0, 0.094d0, 0.2122d0, 0.2729d0, 0.2665d0, 0.3317d0]   
Contains
  Subroutine calGS( )
    Include 'link_fnl_shared.h'
    use LINRG_INT
    Implicit none
    !//ѭ��������itΪ��������ѭ������
    Integer :: it, i
    !//aInit��bInit�ֱ�Ϊ�����²��ʼֵ
    Real(kind=8), parameter :: aInit = 0.1, bInit = 0.1, eps = 1d-8
    !//yestΪ����ֵ��dΪ����ֵ��ʵ��ֵyobs֮��ĲвdTΪ����d�ı��Σ�deltaΪ��������
    Real(kind=8) :: yest(ndata) = 0., d(ndata) = 0., dT(ndata,1) = 0., delta(Nparas,1) = 0.0
    !//JΪ�ſɱȾ���JTΪJ��ת�þ���
    Real(kind=8) :: J(ndata,Nparas) = 0., JT(Nparas,ndata) = 0.
    !//Inv_HΪH�������HΪ��ɭ����
    Real(kind=8) :: H(Nparas,Nparas) = 0., Inv_H(Nparas,Nparas) = 0.
    !//aest��best�ֱ�Ϊ���ݲ���
    Real(kind=8) :: aest = 0.,best = 0.
    Real(kind=8) :: atmp, btmp, tmp
      
    aest = aInit
    best = bInit

    Write(*,"('------------------------------------------------')") 
    Do it = 1, Niters
      yest = aest * xobs / ( best + xobs ) !//���ݵ�ǰaest��best��xobs���õ�����ֵyest
      d = yobs - yest   !//������ֵ֪yobs��yest�Ĳв�
      
      Do i = 1, ndata  !//�����ſɱȾ���dy/da = x / ( b + x )��dy/db = -a*x / ( b + x )**2
        J(i,1) = xobs(i) / ( best + xobs(i) )
        J(i,2) = -aest * xobs(i) / ( best + xobs(i) )**2
      End do
      
      JT = transpose(J)   
      H = matmul( JT,J )  !//���㺣ɭ����
    
      !//���㲽��delta�������ݲ��������µĲ�������ֵ
      call LINRG( H,Inv_H )   !//ʹ��imsl�����⣬����H�������Inv_H
    
      dT = reshape( d,[ndata,1] )  !//Ϊ�������ڲ�����Matmul�ļ��㷨�򣬶�d��������״���иı�
      delta = matmul( matmul( Inv_H,JT ),dT )  !//deltaΪ����
      atmp = aest + delta(1,1)
      btmp = best + delta(2,1)
    
      !//���||delta||<1e-8����ֹ������Ҳ������ǰһ�����һ�ε�aest��best�Ĳ���������
      If ( dot_product(delta(:,1),delta(:,1))<eps ) exit  
      aest = atmp
      best = btmp
    
      Write(*,"('aest =',f10.6,2x,'best =',f10.6)") aest, best
    End Do
    Write(*,"('------------------------------------------------')") 
    Write(*,"('ֹͣ�������ܹ�����',g0,'��')") it - 1
    !//��������������Լ��ٷֱ����������ļ�,�ܹ�100�����ݵ�,��������Ϊ[0-50]
    Open(Fileid,file='fitdat.dat',status='unknown')
    Do i = 0, 100
      tmp = dble(i) / 20.d0
      Write( fileid,'(2f12.8)') tmp, aest * tmp / ( best + tmp )
    End Do
    Close( fileid )
    !//������ݲ���aest��best
    Write(*,"(/,'���ݲ���Ϊ��')")
    Write(*,"('a_est =',f10.6)") aest
    Write(*,"('b_est =',f10.6)") best
  End subroutine calGS
End Module varGS
  

Program GaussNewton
  use varGS
  Implicit None
  call calGS( )
End Program GaussNewton