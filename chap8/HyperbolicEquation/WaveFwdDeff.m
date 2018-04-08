function wave_equation()  %һά������β�������
options={'�ռ�˳�L','�ռ����N' ,'ʱ�����M','�������ٶ�v',...
'�ȶ�������ֵr(ȡֵ����С��1)','��ʼ�ٶȵ�����ʽform(ѡ��1��2)'};
topic='seting';
lines=1;
def={'1','100','2000','1','0.4','1'};
p=inputdlg(options,topic,lines,def);
L=eval(p{1});
N=eval(p{2});
M=eval(p{3});
v=eval(p{4});
r=eval(p{5});%r��ֵ����С��1
form=eval(p{6});
%***************************************************************
h=L/N;%�ռ䲽��
x=0:h:L;
x=x';
tao=r*h/v;%ʱ�䲽��
tm=M*tao;%����������ʱ��tm
t=0:tao:tm;
t=t';
%�����ֵ�ͳ�ֵ
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
%�ò�ַ���Ⲩ������
for j=3:(M+1)
    for i=2:N
        U(i,j)=2*(1-r^2)*U(i,j-1)+r^2*(U(i+1,j-1)+U(i-1,j-1))-U(i,j-2);
    end
end
%���ÿռ�����
for i=1:N+1
    T(i,:)=t;
end
for j=1:M+1
    X(:,j)=x;
end 
%���Ƴ�����ͼ�Σ���U-x-tͼ��
figure(1)
mesh(T,X,U)
xlabel('t');
ylabel('x');
zlabel('U');
title һά�������̵�U-x-tͼ
%���Ƴ�ƽ��ͼ�Σ���U-xͼ��
figure(2)
for k=1:M+1
    plot(x,U(:,k))
    hold on
end 
xlabel('x');
ylabel('U');
title һά�������̵�U-xͼ

function y=border_funo(t)%z=0���ı߽�����
y=0;
return

function y=border_fune(t)%z=L���ı߽�����
y=t.*0;
return

function y=init_fun1(x)%��ʼλ������
%w=3*pi;
%w=pi;
%w=2*pi;
y=0.*x;
return

function y=init_fun2(x)%��ʼ�ٶ�����
y=2*pi*sin(pi*x);
return