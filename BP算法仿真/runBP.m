clc;
clear;
k=1;
img = imread('JellBeans.bmp');
imshow(img)
img = double(img);
subplot(2,3,1);
imagesc(uint8(img)),title('original'); 
xlabel('(a)');
for g=0.4:0.1:0.8
    [mat,p]=CS_BP(g,img);
        subplot(2,3,k+1);
        imagesc(uint8(mat)),title('recovery');
        Str = sprintf('(%s)  M/N=%.1f  %.3fdB',k+97,g,p);
        xlabel(Str);
    k = k+1;
end
