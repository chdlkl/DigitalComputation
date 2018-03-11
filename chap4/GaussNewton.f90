!// 本代码原理以及测试数据来自wiki: https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm#Example
Module varGS
  Implicit None
  !//Nparas为参数个数，ndata为数据个数，Niters为最大迭代次数
  Integer, parameter :: Nparas = 2, ndata = 7, Niters = 50, Fileid = 101
  !//xobs为观测数据自变量，yobs为观测数据因变量，两者均为已知数据
  Real(kind=8), parameter :: xobs(ndata) = [0.038d0, 0.194d0, 0.425d0, 0.626d0, 1.253d0, 2.500d0, 3.740d0]
  Real(kind=8), parameter :: yobs(ndata) = [0.050d0, 0.127d0, 0.094d0, 0.2122d0, 0.2729d0, 0.2665d0, 0.3317d0]   
Contains
  Subroutine calGS( )
    Include 'link_fnl_shared.h'
    use LINRG_INT
    Implicit none
    !//循环变量，it为迭代次数循环变量
    Integer :: it, i
    !//aInit与bInit分别为参数猜测初始值
    Real(kind=8), parameter :: aInit = 0.1, bInit = 0.1, eps = 1d-8
    !//yest为估计值，d为估计值和实际值yobs之间的残差，dT为数组d的变形，delta为增量矩阵
    Real(kind=8) :: yest(ndata) = 0., d(ndata) = 0., dT(ndata,1) = 0., delta(Nparas,1) = 0.0
    !//J为雅可比矩阵，JT为J的转置矩阵
    Real(kind=8) :: J(ndata,Nparas) = 0., JT(Nparas,ndata) = 0.
    !//Inv_H为H的逆矩阵，H为海森矩阵
    Real(kind=8) :: H(Nparas,Nparas) = 0., Inv_H(Nparas,Nparas) = 0.
    !//aest，best分别为反演参数
    Real(kind=8) :: aest = 0.,best = 0.
    Real(kind=8) :: atmp, btmp, tmp
      
    aest = aInit
    best = bInit

    Write(*,"('------------------------------------------------')") 
    Do it = 1, Niters
      yest = aest * xobs / ( best + xobs ) !//根据当前aest，best及xobs，得到函数值yest
      d = yobs - yest   !//计算已知值yobs与yest的残差
      
      Do i = 1, ndata  !//计算雅可比矩阵。dy/da = x / ( b + x )，dy/db = -a*x / ( b + x )**2
        J(i,1) = xobs(i) / ( best + xobs(i) )
        J(i,2) = -aest * xobs(i) / ( best + xobs(i) )**2
      End do
      
      JT = transpose(J)   
      H = matmul( JT,J )  !//计算海森矩阵
    
      !//计算步长delta，并根据步长计算新的参数估计值
      call LINRG( H,Inv_H )   !//使用imsl函数库，计算H的逆矩阵Inv_H
    
      dT = reshape( d,[ndata,1] )  !//为了满足内部函数Matmul的计算法则，对d的数组形状进行改变
      delta = matmul( matmul( Inv_H,JT ),dT )  !//delta为增量
      atmp = aest + delta(1,1)
      btmp = best + delta(2,1)
    
      !//如果||delta||<1e-8，终止迭代。也可以用前一次与后一次的aest与best的差来做条件
      If ( dot_product(delta(:,1),delta(:,1))<eps ) exit  
      aest = atmp
      best = btmp
    
      Write(*,"('aest =',f10.6,2x,'best =',f10.6)") aest, best
    End Do
    Write(*,"('------------------------------------------------')") 
    Write(*,"('停止迭代，总共迭代',g0,'次')") it - 1
    !//输出正反演数据以及百分比误差输出到文件,总共100个数据点,计算区域为[0-50]
    Open(Fileid,file='fitdat.dat',status='unknown')
    Do i = 0, 100
      tmp = dble(i) / 20.d0
      Write( fileid,'(2f12.8)') tmp, aest * tmp / ( best + tmp )
    End Do
    Close( fileid )
    !//输出反演参数aest，best
    Write(*,"(/,'反演参数为：')")
    Write(*,"('a_est =',f10.6)") aest
    Write(*,"('b_est =',f10.6)") best
  End subroutine calGS
End Module varGS
  

Program GaussNewton
  use varGS
  Implicit None
  call calGS( )
End Program GaussNewton