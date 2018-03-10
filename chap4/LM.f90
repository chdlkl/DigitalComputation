!----------------------------------------------------------------------------
!                                     LM算法
!                      说明：程序中的变量和数组在声明时已初始化
!                           方程模型：y = A*Exp( -B*x )
!                           目标：利用LM算法反演参数A与B
!                         原方程：y = 20.5*Exp( -0.24*x )
!                                   作者：luk
!                              完成日期：2017年5月18日
!----------------------------------------------------------------------------
Module VarLM
  Implicit None
  !//Nparas为参数个数，Ndata为数据个数,Niters为最大迭代次数
  Integer     ,Parameter :: Nparas = 2,Ndata = 9,Niters = 50,Fileid = 11
  !//临时变量，order=1，继续迭代，否则停止迭代
  Integer                :: order
  !//循环变量，it为迭代次数循环变量
  Integer                :: it,i,ii
  !//x_1为观测数据自变量，y_1为观测数据因变量，两者均为已知数据
  Real(kind=8),Parameter :: x_1(Ndata) = [0.25,0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0]
  Real(kind=8),Parameter :: y_1(Ndata)  = [19.306,18.182,16.126,14.302,12.685,9.978,7.849,4.857,3.005]   
  !//a0与b0分别为参数猜测初始值,RealA、RealB为真值
  Real(kind=8),Parameter :: a0 = 0.,b0 = 0.,RealA = 20.5,RealB = 0.24
	!//误差限度
  Real(kind=4)           :: eps = 1e-8
  !//y_est为估计值，d为估计值和实际值y_1之间的误差，rd为数组d的变形
  Real(kind=8)           :: y_est(Ndata) = 0.,d(Ndata) = 0.,rd(Ndata,1) = 0.
  !//J为雅可比矩阵，JT为J的转置矩阵，H为海森矩阵
  Real(kind=8)           :: J(Ndata,Nparas) = 0.,JT(Nparas,Ndata) = 0.,H(Nparas,Nparas) = 0.
  !//Inv_H_1m为H_1m的逆矩阵
  Real(kind=8)           :: H_Lm(Nparas,Nparas)=0.,Inv_H_Lm(Nparas,Nparas)=0.
  !//eye为单位矩阵，delta为步长矩阵
  Real(kind=8)           :: eye(Nparas,Nparas) = 0.0,delta(Nparas,1) = 0.0
  !//a_est，b_est分别为反演参数
  Real(kind=8)           :: a_est = 0.,b_est = 0.,e = 0.
  !//计算新对应的值
  Real(kind=8)           :: y_est_Lm(Ndata) = 0.,d_Lm(Ndata) = 0.,e_Lm = 0.
  !//临时变量，lamda为LM算法的阻尼系数初始值
  Real(kind=8)           :: a_Lm,b_Lm,lamda = 0.01,v = 10.d0,temp
 
End Module VarLM
  

Program TestLM
  Include 'link_fnl_shared.h'
  Use VarLM
  Use LINRG_INT
  Implicit None
  
  Write(*,"('------------------------------------------------')") 
  Write(*,"('    方程模型：y = A*Exp( -B*x )     ')")
  Write(*,"('   目标：利用LM算法反演参数A与B     ')")
  Write(*,"('          作者：白南先生            ')") 
  Write(*,"('       完成日期：2017年5月18日      ')") 
  !//生成单位矩阵
  Forall( i = 1:Nparas, ii = 1:Nparas, i==ii ) eye(i,ii) = 1.d0
  
  !//第一步：变量赋值
  order = 1
  a_est = a0
  b_est = b0
  
  !//第二步：迭代
  Write(*,"('------------------------------------------------')") 
  loop1: Do it = 1,Niters
    If ( order==1 ) Then  !//根据当前估计值，计算雅可比矩阵
      
      y_est = a_est*Exp( -b_est*x_1 ) !//根据当前a_est，b_est及x_1，得到函数值y_est
      d = y_1 - y_est   !//计算已知值y_1与y_est的误差
      
      loop2: Do i = 1,Ndata  !//计算雅可比矩阵。dy/da = Exp(-b*x)，dy/db = -a*x*Exp(-b*x)
        J(i,1) = Exp( -b_est*x_1(i) ) 
        J(i,2) = -a_est*x_1(i)*Exp( -b_est*x_1(i) )
      End Do loop2    
      
      JT=Transpose(J)   
      H = Matmul( JT,J )  !//计算海森矩阵
      
    End If
		
    If ( it==1 ) e = Dot_Product(d,d)  !//若是第一次迭代，计算误差epsilon
    H_Lm = H + lamda*eye   !//根据阻尼系数lamda混合得到H矩阵
    
    !//计算步长delta，并根据步长计算新的参数估计值
    Call LINRG( H_Lm,Inv_H_Lm )   !//使用imsl函数库，计算H_Lm的逆矩阵Inv_H_Lm
    
    rd = Reshape( d,[Ndata,1] )  !//为了满足内部函数Matmul的计算法则，对d的数组形状进行改变
    delta = Matmul( Inv_H_Lm,matmul( JT,rd ) )  !//delta为增量
    a_Lm = a_est + delta(1,1)
    b_Lm = b_est + delta(2,1)
    
    !//如果||delta||<1e-8，终止迭代
    If ( Dot_Product(delta(:,1),delta(:,1))<eps ) Exit
		
    !//计算新的可能估计值对应的y和计算残差e
    y_est_Lm = a_Lm*Exp( -b_Lm*x_1 )
    d_Lm = y_1 - y_est_Lm
    e_Lm = Dot_Product( d_Lm,d_Lm )  !//e_LM等于||y_1 - y_est_LM||
    
    !//根据误差，决定如何更新参数和阻尼系数
    !//迭代成功时将lamda减小，否则增大lamda
    If ( e_Lm<e ) Then
      lamda = lamda/v
      a_est = a_Lm
      b_est = b_Lm
      e = e_Lm
      order = 1
    Else
      order = 0
      lamda = lamda*v
    End If
    
    Write(*,"('a_est=',g0,2X,'b_est=',g0)") a_est,b_est
  End Do loop1
  Write(*,"('------------------------------------------------')") 
  Write(*,"('停止迭代，总共迭代',g0,'次')") it - 1
  !//输出正反演数据以及百分比误差输出到文件,总共100个数据点,计算区域为[0-50]
  Open(Fileid,file='原始数据与反演数据.dat',status='unknown')
  Do i = 0,100
    temp = i/2.d0
    Write(Fileid,'(f7.3,3f15.8)') temp,RealA*Exp( -RealB*temp ),a_est*Exp( -b_est*temp ),&
      Abs( a_est*Exp( -b_est*temp )-RealA*Exp( -RealB*temp ) )/RealA/Exp( -RealB*ii )*100
  End Do
  Close(Fileid)
  !//输出反演参数a_est，b_est以及原始数据和拟合数据
  Write(*,"(/,'反演参数为：')")
  Write(*,"('a_est=',g0)") a_est
  Write(*,"('b_est=',g0)") b_est
  Write(*,"(/,'原始数据为：')") 
  Write(*,'(3f9.4/)') y_1
  Write(*,"('拟合数据为：')") 
  Write(*,'(3f9.4/)') a_est*Exp( -b_est*x_1 )
  Write(*,"('------------------------------------------------')") 
  
End Program TestLM