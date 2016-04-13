function [img_mat,PSNR] = CS_BP(g,picture)

% ����ͼ��
% picture=imread('Cameraman.bmp');    
% imshow(picture);
img=double(picture);
[img_height,img_width]=size(img);


% ���ɸ�˹�����������
Phi=Get_Gauss_Mat(floor(img_height/(1/g)),img_width); 


% ���DCT�任����
Dct_Mat=dctmtx(256)';

% Y = ��* X, ��ò�������Y
Y_Mat=Phi*img;        

% ����BP�㷨�ָ�X
Sparse_Mat=zeros(img_height,img_width);  
Cmpressed_Mat=Phi*Dct_Mat;
for i=1:img_width
    column_rec=bp(Y_Mat(:,i),Cmpressed_Mat);
    Sparse_Mat(:,i)=column_rec';           
end
img_recovery=Dct_Mat*Sparse_Mat;         
psnr = 20*log10(255/sqrt(mean((img(:)-img_recovery(:)).^2)));
PSNR = psnr;
img_mat = img_recovery;

% ���չʾ

% subplot(1,2,1);
% imagesc(uint8(img)),title('original');
% xlabel('(a)');
% subplot(1,2,2);
% imagesc(uint8(img_recovery)),title('recovery')
% xlabel(strcat('(b)',' ',num2str(psnr),'dB'));
% figure(2);
% subplot(1,2,1),imagesc(Phi),title('measurement mat')
% subplot(1,2,2),imagesc(Dct_Mat),title('dct mat')

