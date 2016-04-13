clc;
clear;
k=1;
img = imread('lena.bmp');
img = double(img);
psnr_BP=zeros(1,8);
for g=0.1:0.1:0.8
    [mat,p]=CS_BP(g,img);
        psnr_BP(k)=p;
    k = k+1;
end
