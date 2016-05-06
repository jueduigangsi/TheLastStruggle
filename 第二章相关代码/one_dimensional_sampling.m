clc;
clear;
function ONE_DIM_SAMP()
    k = 1:100;
    N = 100;
    y = sin(4*k*pi/N);
    y1 = dct2(y);
    figure;
    stem(y,'.');
    title('y = sin(4*k*pi/N),k=1:100,N=100');
    figure;
    stem(y1,'.');
    title('DCT( y )');
end