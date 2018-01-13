!// 四类三次样条插值
!// 1.自然三次样条插值; 2.钳制三次样条插值; 3.曲率(二阶导数)任意调整的三次样条插值; 4.抛物线端点的三次样条曲线; 5.非纽结三次样条(matlab中的命令为spline)
Module CalSpline
  Implicit none 
  Character(*), parameter :: file1 = 'xy.txt'  !// 插值基点文件
  Real(kind=8), allocatable :: x(:), y(:)  !// 已知节点
  Real(kind=8), allocatable :: coeff(:,:)  !// 系数矩阵b,c,d
  Real(kind=8), allocatable :: a(:,:), r(:,:) !// a为求解系数c的左端项系数矩阵, r为求解系数c的右端项
  Real(kind=8), allocatable :: segma(:), delta(:)  !// segma(i) = x(i+1) - x(i), delta(i) = y(i+1) - y(i), i = 1,...,numLine-1
  Integer :: numLine
Contains
Subroutine CalFileLine ()
  Implicit none
  Integer, parameter :: fileid = 101
  Integer :: info = 0, i
  !//---------读取节点数据个数---------
  numLine = 0
  Open ( fileid, file = file1 )
  Do 
    Read ( fileid, fmt = *, iostat = info )
    If ( info /= 0 ) Exit
    numLine = numLine + 1
  End Do 
  print*, 'numLine:'
  print*, numLine
  Allocate( x(numLine), y(numLine), stat = info )
  If ( info == 0 ) Then 
    Write ( *,'(1x,g0)' ) "Allocate x, y array successfully!"
  Else
    Write ( *,'(1x,g0)' ) "Allocate x, y array fail!"
  End if 
  !//---------读取xy---------
  rewind ( fileid )
  print*, 'array xy:'
  Do i = 1, numLine
    Read ( fileid,* ) x(i), y(i)
    print*, x(i), y(i)
  End Do 
  Close ( fileid )

End subroutine CalFileLine

Subroutine Calcoeff ()
  Use lapack95  !// 使用lapack函数库求解逆矩阵
  Implicit none 
  Integer :: info = 0, i
  Real(kind=8) :: tmp_a(numLine,numLine), tmp_r(numLine,1)
  Real(kind=8) :: v1 = 0.d0, vn = 0.d0  !// v1, vn是两个端点的一阶导数，由用户自行设定，用于钳制三次样条插值.(这里为了方便，v1,vn在代码里也可以用于曲率调整的三次样条)
  !//---------分配数组---------
  Allocate( segma(numLine-1), delta(numLine-1), coeff(numLine,3), a(numLine,numLine), r(numLine,1), stat = info )
  segma = 0.d0; delta = 0.d0
  coeff = 0.d0; a = 0.d0; r = 0.d0
  If ( info == 0 ) Then 
    Write ( *,'(1x,g0)' ) "Allocate array successfully!"
  Else
    Write ( *,'(1x,g0)' ) "Allocate array fail!"
  End if 
  !//---------计算segma与delta---------
  Do i = 1, numLine - 1
    segma(i) = x(i+1) - x(i)
    delta(i) = y(i+1) - y(i)
  End Do 
  !//---------构建求取系数c的左端项矩阵a与右端项r---------
  ! a(1,1) = 1.d0; a(numLine,numLine) = 1.d0
  ! r(1,1) = 0.d0; r(numLine,1) = 0.d0  !// 自然三次样条Nature cubic spline: 两端点处的二阶导数为0
  ! a(1,1:2) = [ 2.d0 * segma(1), segma(1) ]; a(numLine,numLine-1:numLine) = [ segma(numLine-1), 2.d0 * segma(numLine-1) ]
  ! r(1,1) = 3.d0 * ( delta(1) / segma(1) - v1 ); r(numLine,1) = 3.d0 * ( vn - delta(numLine-1) / segma(numLine-1) )  !// 钳制三次样条插值，两端点处的一阶导数为0
  ! a(1,1) = 2.d0; a(numLine,numLine) = 2.d0
  ! r(1,1) = v1; r(numLine,1) = vn  !// 曲率(二阶导数)任意调整的三次样条插值，这里的v1,vn为两端点的二阶导数(曲率)
  ! a(1,1:2) = [ 1.d0, -1.d0 ]; a(numLine,numLine-1:numLine) = [ 1.d0, -1.d0 ]
  ! r(1,1) = 0.d0; r(numLine,1) = 0.d0  !// 抛物线端点的三次样条曲线，通过使三次项的系数为0，得样条的起始和结束部分S1和Sn-1至多2阶
  a(1,1:3) = [ segma(2), -( segma(1) + segma(2) ), segma(1) ]; a(numLine,numLine-2:numLine) = [ segma(numLine-1), -( segma(numLine-2) + segma(numLine-1) ), segma(numLine-2) ]
  r(1,1) = 0.d0; r(numLine,1) = 0.d0  !// 非纽结三次样条(X2,Xn-1为非纽结点)
  Do i = 2, numLine - 1
    a( i,i-1:i+1 ) = [ segma(i-1), 2.d0* ( segma(i-1) + segma(i) ), segma(i) ]  !// 左端项系数矩阵
    r( i,1 ) = 3.d0 * ( delta(i) / segma(i) - delta(i-1) / segma(i-1) )  !// 右端项系数矩阵
  End Do
  !//---------计算系数c---------
  tmp_a = a; tmp_r = r 
  call gesv( tmp_a, tmp_r )  !// 使用lapack函数库
  coeff(:,2) = tmp_r(:,1)
  !//---------计算系数b,d---------
  Do i = 1, numLine - 1
    coeff(i,1) = delta(i) / segma(i) - segma(i) * ( 2.d0*coeff(i,2) + coeff(i+1,2) ) / 3.d0
    coeff(i,3) = ( coeff(i+1,2) - coeff(i,2) ) / 3.d0 / segma(i)
  End Do 
  print*, 'coeff_b:'
  print*, coeff(:,1)
  print*, 'coeff_c:'
  print*, coeff(:,2)
  print*, 'coeff_d:'
  print*, coeff(:,3)

End subroutine Calcoeff

Subroutine  CalInterpolation ()
  Implicit none 
  Character(*), parameter :: filename = 'xx.txt'  !// 待插值节点文件
  Integer :: info = 0, n, i, j
  Real(kind=8) :: dx
  Real(kind=8), allocatable :: xx(:), yy(:)
  
  !//---------读取待插值节点---------
  n = 0
  Open ( 101, file = filename )
  Do
    Read ( 101, fmt = *, iostat = info ) 
    If ( info /= 0 ) Exit
    n = n + 1
  End Do
  Allocate( xx(n), yy(n), stat = info )
  If ( info == 0 ) Then 
    Write ( *,'(1x,g0)' ) "Allocate xx, yy array successfully!"
  Else
    Write ( *,'(1x,g0)' ) "Allocate xx, yy array fail!"
  End if 
  rewind( 101 )
  print*, 'n:'
  print*, n
  Read ( 101,* ) xx
  Close ( 101 )
  
  !//---------求取插值---------
  outdo:Do i = 1, numLine - 1
    indo:Do j = 1, n 
      If ( xx(j) >= x(i) .and. xx(j) <= x(i+1) ) Then
        dx = xx(j) - x(i)
        yy(j) = coeff(i,3) * dx  !// 使用嵌套乘法求值
        yy(j) = ( yy(j) + coeff(i,2) ) * dx
        yy(j) = ( yy(j) + coeff(i,1) ) * dx + y(i)
      End if 
    End Do indo
  End Do outdo
  !//---------输出插值节点数据---------
  Open ( 101, file = 'Interpolation5.dat' )
  Do i = 1, n 
    Write ( 101,'(*(2x,g0))' ) xx(i), yy(i)
  End Do 
  Close ( 101 )
  Deallocate( xx, yy )
End subroutine CalInterpolation 

Subroutine DeallocateArry ()
  Implicit none
  Integer :: info = 0
  !//---------释放内存---------
  Deallocate( x, y, coeff, a, r, segma, delta, stat = info )
  If ( info == 0 ) Then 
    Write ( *,'(1x,g0)' ) "Deallocate all array successfully!"
  Else
    Write ( *,'(1x,g0)' ) "Deallocate all array fail!"
  End if 
End subroutine DeallocateArry

End module CalSpline

Program spline 
  use CalSpline
  Implicit none
  Call CalFileLine ()
  Call Calcoeff ()
  Call CalInterpolation ()
  Call DeallocateArry ()
End program spline  
