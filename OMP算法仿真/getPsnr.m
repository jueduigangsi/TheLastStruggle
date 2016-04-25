
k=1;
img = imread('lena.bmp');
img = double(img);
% psnr_OMP=zeros(1,8);
for g=4:256
    [mat,p]=CS_OMP(g,img);
        psnr_OMP(g)=p;
    k = k+1;
end
