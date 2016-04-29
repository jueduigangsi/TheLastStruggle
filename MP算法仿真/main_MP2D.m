clc
close all
clear all
I = imread('lena128.bmp');
[m,n] = size(I);
N = m * n;  %稀疏信号长度
vectorI = reshape(I,N,1);
x0 = dct2(vectorI);
sampleRate = 0.3;  %采样率 = M/N
M = round( N * sampleRate );
%------  高斯感知矩阵Phi   -------------
Phi = sqrt(1/M) * randn(M,N);
for i = 1:N
    Phi(:,i) = Phi(:,i) / norm(Phi(:,i));
end
%-------- 测量向量 y  ----------
y  = zeros(M,1);
y = Phi * x0;

%% -----2. MP Reconstruction ------------
x = zeros(N,1);       %x0的逼近信号x
times =  M/2        %迭代次数 = 稀疏度
g = zeros(N,1);       %余量和感知矩阵内积
r  = y;               %余量初始化为y
for count=0:times
    count;
    g = Phi' * r;
    [val,K] = max( abs(g) ) ;
    x(K,1) = x(K,1) + g(K,1);
    r = r - g(K,1) * Phi(:,K);
end
idct_x = idct(x);
rec_x = reshape(idct_x,m,n);
%% --------3.重建质量评价--------------------
relative_err = norm(x-x0) / norm(x)
MSE = norm(x-x0) / (m*n)
SNR = 10*log10( norm(x)^2 / MSE )
Mat_rate = 1 - norm(abs(x)-abs(x0))/norm(abs(x)+abs(x0))
figure
subplot(121),imshow(I,[]),title('orignal image');
subplot(122),imshow(rec_x,[]), title('reconstruct image');


%%  ----测试数据记录--------------------------
% times = M,sampleRate = 0.5;
% relative_err =
%     0.1278
% MSE =
%     0.2658
% SNR =
%    84.3585

% times = M/2,sampleRate = 0.3;
% 迭代次数为614，重建信号稀疏度K=588
% relative_err =
%     0.1593
% MSE =
%     0.3299
% SNR =
%    83.3882
% Mat_rate =
%     0.9235
