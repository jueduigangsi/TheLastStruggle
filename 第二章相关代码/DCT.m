A=imread('lena.bmp');
C=dct2(A); %进行余弦变换
figure;
B=log(abs(C));

figure;
imshow(A);
title('原始图像');
xlabel('(a)');

C(abs(C)<10)=0; %将DCT变换后的系数值小于10的元素设为0
%E=idct2(C);
D=idct2(C)./255; %对DCT变换值归一化，进行余弦反变换???
figure;
imshow(D) ;
title('重构图像');
xlabel('(b)');
% imshow(uint8(E)); is the same as D=idct2(C)./255
% imshow(E,[]); is the same as D=idct2(C)./255

FF=abs(C)<10; %Compute the number of elements which are smaller than 10
sum(sum(FF)) %result:46103
GG=abs(C)>10; %Compute the number of elements which are larger than 10
sum(sum(GG)) %result:19433
