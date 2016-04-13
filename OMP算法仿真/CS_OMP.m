function [img_mat,PSNR] = CS_OMP(g,picture)

% img=imread('lena.bmp');   
% imshow(img)
img=double(picture);
[img_height,img_width]=size(img);

% ���ɸ�˹�����������
Phi=Get_Gauss_Mat(floor(img_height/(1/g)),img_width); 

% ���DCT�任����
Dct_Mat=dctmtx(256)';


% Y = ��* X, ��ò�������Y
Y_Mat=Phi*img;  


% ����SP�㷨�ָ�X
Sparse_Mat=zeros(img_height,img_width);  
Cmpressed_Mat=Phi*Dct_Mat;
for i=1:img_width
    column_rec=omp(Y_Mat(:,i),Cmpressed_Mat,img_height);
    Sparse_Mat(:,i)=column_rec';           
end
img_recovery=Dct_Mat*Sparse_Mat;         
psnr = 20*log10(255/sqrt(mean((img(:)-img_recovery(:)).^2)));
PSNR = psnr;
img_mat = img_recovery;

% ���չʾ
% figure(1)
% subplot(2,2,1),imagesc(img),title('original image')
% subplot(2,2,2),imagesc(Phi),title('measurement mat')
% subplot(2,2,3),imagesc(mat_dct_1d),title('1d dct mat')
% psnr = 20*log10(255/sqrt(mean((img(:)-img_rec_1d(:)).^2)))
% subplot(2,2,4),imagesc(img_rec_1d),title(strcat('1d rec img ',num2str(psnr),'dB'))
% 
% disp('over')








