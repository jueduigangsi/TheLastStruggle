A=imread('lena.bmp');
C=dct2(A); %�������ұ任
figure;
B=log(abs(C));

figure;
imshow(A);
title('ԭʼͼ��');
xlabel('(a)');

C(abs(C)<10)=0; %��DCT�任���ϵ��ֵС��10��Ԫ����Ϊ0
%E=idct2(C);
D=idct2(C)./255; %��DCT�任ֵ��һ�����������ҷ��任???
figure;
imshow(D) ;
title('�ع�ͼ��');
xlabel('(b)');
% imshow(uint8(E)); is the same as D=idct2(C)./255
% imshow(E,[]); is the same as D=idct2(C)./255

FF=abs(C)<10; %Compute the number of elements which are smaller than 10
sum(sum(FF)) %result:46103
GG=abs(C)>10; %Compute the number of elements which are larger than 10
sum(sum(GG)) %result:19433
