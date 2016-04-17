clc;
clear;
k = 1:200;
N = 200;
y = sin(8*k*pi/N);
stem(y,'.');
title('y = sin(8*k*pi/N),k=1:200,N=200');