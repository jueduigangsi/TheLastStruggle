%获取DCT变换基矩阵
function Dct_Mat = Get_Dct_Basis()
    Dct_Mat = zeros(256,256);  
    for m=0:255 
        dcttmp=cos([0:255]'*m*pi/256);
        if m>0
            dcttmp=dcttmp-mean(dcttmp); 
        end;
        Dct_Mat(:,m+1)=dcttmp/norm(dcttmp);
    end