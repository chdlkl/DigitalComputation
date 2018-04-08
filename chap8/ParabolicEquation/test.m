% a = load('output1.dat');
% x = (0:10)*0.1;
% t = (0:250)*0.004;
% w = a(:,3);
% w = reshape( w, 11, 251 );
% mesh( x, t, w' );
% view( 60,30 )
% axis( [0 1 0 1 -1 1] )
function w = test( xl, xr, yb, yt, M, N )
f = @(x) sin(2*pi*x).^2;
l = @(t) 0*t;
r = @(t) 0*t;
D = 1;
h = (xr-xl)/M; k = (yt-yb)/N; m = M - 1; n = N;
sigma = D*k/h/h;
a = diag( 1-2*sigma*ones(m,1) ) + diag( sigma*ones(m-1,1),1 );
a = a + diag( sigma*ones(m-1,1),-1 );
lside = l( yb+(0:n)*k );
rside = r( yb+(0:n)*k );
w(:,1) = f( xl+(1:m)*h )';
for j = 1: n
    w(:,j+1) = a*w(:,j) + sigma*[lside(j);zeros(m-2,1);rside(j)];
end
w = [lside;w;rside]
x = (0:m+1)*h; t = (0:n)*k;
mesh( x, t, w' );
view( 60,30 )
axis( [0 1 0 1 -1 1] )
