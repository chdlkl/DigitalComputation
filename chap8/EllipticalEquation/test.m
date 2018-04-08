a = load('output.dat');
x = (0:0.1:1);
t = (0:0.1:1);
w = a(:,4);
w = reshape( w, 11, 11 );
mesh( x, t, w' );
view( 60,30 )
axis( [0 1 0 1 0 1] )