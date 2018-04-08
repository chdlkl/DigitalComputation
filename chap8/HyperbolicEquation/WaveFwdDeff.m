function wave_equation()  %一维线性齐次波动方程
options={'空间杆长L','空间点数N' ,'时间点数M','波的相速度v',...
'稳定条件的值r(取值必须小于1)','初始速度调用形式form(选择1或2)'};
topic='seting';
lines=1;
def={'1','100','2000','1','0.4','1'};
p=inputdlg(options,topic,lines,def);
L=eval(p{1});
N=eval(p{2});
M=eval(p{3});
v=eval(p{4});
r=eval(p{5});%r的值必须小于1
form=eval(p{6});
%***************************************************************
h=L/N;%空间步长
x=0:h:L;
x=x';
tao=r*h/v;%时间步长
tm=M*tao;%波传播的总时间tm
t=0:tao:tm;
t=t';
%计算边值和初值
U=zeros(N+1,M+1);
Uo=border_funo(t);
Ue=border_fune(t);
Ui=init_fun1(x);
dUi=init_fun2(x);
U(1,:)=Uo;
U(N+1,:)=Ue;
U(:,1)=Ui;
if form==1
  U(:,2)=init_fun1(x)+tao*init_fun2(x);
else
  for i=2:N
  U(i,2)=(1-r^2)*Ui(i)+0.5*r^2*(Ui(i+1)+Ui(i-1))+tao*dUi(i);
  end
end    
%用差分法求解波动方程
for j=3:(M+1)
    for i=2:N
        U(i,j)=2*(1-r^2)*U(i,j-1)+r^2*(U(i+1,j-1)+U(i-1,j-1))-U(i,j-2);
    end
end
%设置空间网格
for i=1:N+1
    T(i,:)=t;
end
for j=1:M+1
    X(:,j)=x;
end 
%绘制出立体图形，即U-x-t图形
figure(1)
mesh(T,X,U)
xlabel('t');
ylabel('x');
zlabel('U');
title 一维波动方程的U-x-t图
%绘制出平面图形，即U-x图形
figure(2)
for k=1:M+1
    plot(x,U(:,k))
    hold on
end 
xlabel('x');
ylabel('U');
title 一维波动方程的U-x图

function y=border_funo(t)%z=0处的边界条件
y=0;
return

function y=border_fune(t)%z=L处的边界条件
y=t.*0;
return

function y=init_fun1(x)%初始位移条件
%w=3*pi;
%w=pi;
%w=2*pi;
y=0.*x;
return

function y=init_fun2(x)%初始速度条件
y=2*pi*sin(pi*x);
return