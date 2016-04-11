A=imread('lena.bmp');
imshow(A) %A is unit8(0,255)
C=dct2(A); %�������ұ任
figure;
B=log(abs(C));
imshow(B)
colormap(jet(64)); %��ʾΪ64���Ҷ�
colorbar; %��ʾ��ɫ������ʾ�任���ϵ���ֲ�
C(abs(C)<20)=0; %��DCT�任���ϵ��ֵС��10��Ԫ����Ϊ0
%E=idct2(C);
D=idct2(C)./255; %��DCT�任ֵ��һ�����������ҷ��任???
figure;
imshow(D) ;
% imshow(uint8(E)); is the same as D=idct2(C)./255
% imshow(E,[]); is the same as D=idct2(C)./255

FF=abs(C)<10; %Compute the number of elements which are smaller than 10
sum(sum(FF)) %result:56632
GG=abs(C)>10; %Compute the number of elements which are larger than 10
sum(sum(GG)) %result:16025
