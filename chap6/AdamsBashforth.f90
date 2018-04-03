!---------------------------------------------------------
!             此程序用的是AdamsBashforth(多步外插法)求解
!                       y'=-3*y
!                       y(0) = 1
!                       t = [0,2]
!                   解析解：y = e^(-3t)
!---------------------------------------------------------
Module AdamsBashforth
  Implicit none
Contains
  Subroutine solve()
    Implicit none
    Real(kind=8), parameter :: h = 1d-2, T0 = 0.d0, T1 = 2.d0
    Integer, parameter :: n = 3, num = nint( (T1-T0)/h )  !// n=3,(0,1,2,3)表示四阶四步法
    Real(kind=8) :: a(0:n), b(0:n), sum, t, sum1, b1, jL, u(0:n), u1(0:num)
    Integer :: i, j, L
  
    Write(*,*) "----------求解Adams外插公式系数值-----------"
    a(0) = 1.d0
    Do i = 1, n
      sum = 0.d0; t = 2.d0
      Do j = i, 1, -1
        sum = sum + a(j-1) / (t)
        t = t + 1.d0
      End do
      a(i) = 1.d0 - sum
    End do
    Write(*,*) a
  
    Write(*,*) "----------验证Adams外插公式系数值-----------"
    Print*
    !//验证Adams外插公式系数值分解是否正确
    sum1 = 0.d0
    Do i = 0, n
      sum1 = sum1 + a(n-i) / (1+i)
    End do
    If( abs(sum1-1.0) < 1d-6 ) then   !//sum1=1.0,说明Adams外插公式系数分解正确
      write(*,*) "----------Adams外插公式系数分解正确---------"  
    else
      write(*,*) "---------Adams外插公式系数分解不正确--------"
      stop
    End If
  
    Print*, '--------------------------------------------'
    Write(*,*) "---------------求解系数Bkl值----------------"
    Print*, '--------------------------------------------'
    Write(*,*) 
  
    Do L = 0, n
      b1 = 0.0
      Do j = L, n
        call calculate_jL(j,L,jL)  !//求解(s(s-1)...(s-j+1)/j!)
        b1 = b1 + (-1)**L * a(j) * jL
      End do
      b(L) = b1
    End do
    Write(*,*) b*24
  
    !//由梯形算法计算前四个解
    call Trapezoid( n, u1, num, h, t0 )
  
    !//AdamsForth公式
    Do i = n+1, num
      u1(i) = u1(i-1) + h * ( b(0)*(-3.*u1(i-1)) + b(1)*(-3.*u1(i-2)) + b(2)*(-3.*u1(i-3)) + b(3)*(-3.*u1(i-4)) )
    End do
  
    open( 100, file = 'waicha.dat', status = 'unknown' )
    Do i = 0, num
      write( 100,'(4f12.8)' ) i*h, u1(i), exp(-3.d0*h*i), abs((u1(i)-exp(-3.d0*h*i)) / exp(-3.d0*h*i)) * 100
    End do
    close( 100 ) 
  End subroutine solve
  
  Subroutine calculate_jL( j, L, jL )
    Implicit none
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
        Product1 = Product1 * dble(i)  !//求解分子s(s-1)...(s-j+1)
      End Do
      Do i = j, j-L+1, -1
        Product2 = Product2 * dble(i)  !//求解j!
      End Do
      jL = Product2 / Product1
    End If
  End Subroutine calculate_jL
  
  Subroutine Trapezoid( n, u1, num, h, t0 )
    Implicit none
    Integer, intent(in) :: n, num
    Real(kind=8), intent(inout) :: u1(0:n)
    Real(kind=8), intent(in) :: h, t0
    Real(kind=8) :: y0, t, y
    Integer :: i
    
    y0 = 1.d0  !// 初值
    u1(0) = y0
    Do i = 1, n
      t = t0 + h * ( i - 1 )
      !// Wi+1 = Wi + h * ( f(Ti,Wi) + f( Ti+h,Wi+h*f(Ti,Wi) ) )
      y = y0 + h * ( derfunc( y0 ) + derfunc( y0 + h*derfunc( y0 ) ) ) / 2.d0
      y0 = y
      u1(i) = y
    End do
  End subroutine Trapezoid
  
  Real(kind=8) function derfunc( y )
    Implicit none
    Real(kind=8), intent(in) :: y
    derfunc = -3.d0 * y  !// 已知的微分方程
  End function derfunc
  
End module AdamsBashforth
  
  
Program main
  use AdamsBashforth
  Implicit none
  call solve()
End Program main