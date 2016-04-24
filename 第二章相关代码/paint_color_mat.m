clear;clc;
A=eye(12);
p=randperm(100);

 B=A;B(~A)=nan;
pcolor(B);colorbar;
axis ij;