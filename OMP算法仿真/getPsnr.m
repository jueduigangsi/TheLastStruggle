clc;
clear;
k=1;
img = imread('lena.bmp');
img = double(img);
psnr_OMP=zeros(1,8);
for g=0.1:0.1:0.8
    [mat,p]=CS_OMP(g,img);
        psnr_OMP(k)=p;
    k = k+1;
end
