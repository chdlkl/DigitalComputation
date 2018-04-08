a = load('output.dat');
x = (0:100)*0.01;
t = (0:2000)*(1/2000);
w = a(:,3);
w = reshape( w, 101, 2001 );
mesh( x, t, w' );
view( 60,30 )
% axis( [0 1 0 1 -1 1] )