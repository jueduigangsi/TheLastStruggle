function X=omp(Y,Compressed_Mat,m)

n=length(Y);
s=floor(n/4);                                     %  测量值维数
X=zeros(1,m);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=Y;                                            %  残差值 

for times=1:s;                                  %  迭代次数(稀疏度是测量的1/4)

    product=abs(Compressed_Mat'*r_n);
    
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,Compressed_Mat(:,pos)];                   %  矩阵扩充
    Compressed_Mat(:,pos)=zeros(n,1);                      %  选中的列置零（实质上应该去掉，为了简单将其置零）
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*Y;           %  最小二乘,使残差最小
    r_n=Y-Aug_t*aug_x;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    
end
X(pos_array)=aug_x;                           %  重构的向量 


