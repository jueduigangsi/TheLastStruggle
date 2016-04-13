clc;
clear;
k=1;
img = imread('lena.bmp');
img = double(img);
psnr_SP=zeros(1,8);
for g=0.1:0.1:0.8
    [mat,p]=CS_SP(g,img);
        psnr_SP(k)=p;
    k = k+1;
end
