!---------------------------------------------------------
!             此程序用的是AdamsBashforth(多步外插法)求解
!                       y'=-3*y
!                       y(0) = 1
!                       t = [0,2]
!                   解析解：y = e^(-3t)
!---------------------------------------------------------
Module AdamsMoulton
  Implicit none
Contains
  Subroutine solve()  
    Real(kind=8), parameter :: T0 = 0.d0, T1 = 2.d0 !//n为阶数，T1，T0为计算区间
    Real(kind=8), parameter :: h = 1d-2  !//h为步长
    Integer, parameter :: n = 3
    Integer, parameter :: num = nint((T1-T0)/h)  !//num为总点数
    Real(kind=8) :: a(0:n), b(0:n), sum, t, sum1, b1, jL, u1(0:num), t11, t22, u11(0:num), tmp
    !//数组a存放生成函数法的导出系数，数组b存放最终表达式前的系数
    !//sum,sum1,b1都为初始变量，用来存放循环值，t11,t22为临时变量，进行迭代收敛
    !//u1存放最终数值解，u11存放用外插法计算的数值解
    Integer :: i, j, L
  
    write(*,*) "----------求解Adams内插公式系数-------------"
    a(0) = 1.d0
    Do i = 1, n
      sum = 0.d0; t = 2.d0
      Do j = i, 1, -1
        sum = sum + a(j-1) / (t)
        t = t + 1.d0
      End Do
      a(i) = 0.d0 - sum  !//此处与外插法不一样，a(i)=0.-sum;而外插法为a(i)=1.-sum
    End Do
    write(*,*) a
  
    Print*,"----------验证Adams内插公式系数-------------"
    Print*
    sum1 = 0.d0  !//验证Adams内插公式系数值分解是否正确
    Do i = 0, n
      sum1 = sum1 + a(n-i) / dble(1+i)
    End do
    If( abs(sum1-0.0) < 1e-6 ) then   !//sum1=0.0,说明Adams内插公式系数分解正确
      write(*,*) "----------Adams内插公式系数分解正确---------"  
    else
      write(*,*) "---------Adams外插公式系数分解不正确--------"
      stop
    End if
  
    Print*, '--------------------------------------------'
    Write(*,*) "---------------求解系数Bkl值----------------"
    Print*, '--------------------------------------------'
  
    Do L = 0, n
      b1 = 0.d0
      Do j = L, n
        call calculate_jL( j, L, jL )  !//求解(s(s-1)...(s-j+1)/j!)
        b1 = b1 + (-1)**L * a(j) * jL
      End do
      b(L) = b1
    End do
    write(*,*) b*24  !//24是为了将系数显示成整数，没有其他特别含义  

    open( 100, file = 'waicha.dat', status = 'old')  !//读取Adams外插法计算的数值
    Do i = 0, num
      read(100,*) tmp, u11(i)
      If( i<n ) u1(i) = u11(i)   !//前3个数值由外插算法给出
    End do
    close( 100 )
  
    !//AdamsMoulton公式
    Do i = n, num  !由于本程序用的是3阶，所以内插法的前3个数值要用外插法给出（内插法要给出前4个数值）
      t11 = u11(i)
      Do
        t22 = u1(i-1) + h * ( b(0)*(-3.*t11) + b(1)*(-3.*u1(i-1)) + b(2)*(-3.*u1(i-2)) + b(3)*(-3.*u1(i-3)) )
        If( abs(t22-t11) < 1d-6 ) then  
          exit  !//如果迭代收敛的话，跳出循环，赋值给u1(i)
        else
          t11 = t22  !//不满足迭代要求，继续迭代
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
        Product1 = Product1 * i  !//求解分子s(s-1)...(s-j+1)
      End do
      Do i = j, j-L+1, -1
        Product2 = Product2 * i  !//求解j的阶乘
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
  
