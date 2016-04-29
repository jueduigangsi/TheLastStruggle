clc;
clear;
k=1;
img_str={'lena.bmp','Clock.tiff','JellBeans.bmp','Moonsurface.tiff','house.bmp','Cameraman.bmp'};

g=0.8;
for k=1:3
    img = imread(char(img_str(k)));
    if k==1
      imshow(img)
    end
    img = double(img);
    subplot(2,3,k);
    imagesc(uint8(img)),title('original'); 
    Str1 = sprintf('(%s)',k+96);
    xlabel(Str1);

        [mat,p]=CS_OMP(g,img);
        subplot(2,3,k+3);
        imagesc(uint8(mat)),title('recovery');
        Str = sprintf('(%s)  M/N=%.1f  %.3fdB',k+96+3,g,p);
        xlabel(Str);
    
end
