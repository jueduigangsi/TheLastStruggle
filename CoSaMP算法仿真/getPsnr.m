clc;
clear;
k=1;
img = imread('lena.bmp');
img = double(img);
psnr_CoSaMP=zeros(1,8);
for g=0.1:0.1:0.8
    [mat,p]=CS_CoSaMP(g,img);
        psnr_CoSaMP(k)=p;
    k = k+1;
end
