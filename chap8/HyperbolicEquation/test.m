a = load('output.dat');
x = (0:0.1:1);
t = (1:0.1:2);
w = a(:,3);
w = reshape( w, 11, 11 );
mesh( x, t, w' );
view( 60,30 )
axis( [0 1 1 2 0 2] )