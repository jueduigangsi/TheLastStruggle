function DCT()
    img=imread('lena.bmp');
    img_data=dct2(img);
    figure;
    imshow(img);
    title('original');
    xlabel('(a)');
    img_data(abs(img_data)<10)=0; 
    img_rec=idct2(img_data)./255; 
    figure;
    imshow(img_rec) ;
    title('recovery');
    xlabel('(b)');
end