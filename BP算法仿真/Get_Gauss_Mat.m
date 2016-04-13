function Gauss_Mat = Get_Gauss_Mat( M, N )
%GET_ GAUSS_MAT 生成随机高斯矩阵
Gauss_Mat=randn(M,N); 
Gauss_Mat = Gauss_Mat./repmat(sqrt(sum(Gauss_Mat.^2,1)),[M,1]);



