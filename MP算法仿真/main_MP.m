clc
clear all;
close all;
M = 800;   %观测信号长度
N = 256*256;  %稀疏信号长度
K = 23;   %稀疏度
%% -----1.生成稀疏度为K的稀疏信号-----
x0 = zeros(N,1);
p = randperm(N);
x0(p(1:K),1) = randn(K,1);
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
times =  K;           %迭代次数 = 稀疏度
g = zeros(N,1);       %余量和感知矩阵内积
r  = y;               %余量初始化为y
for n=0:times
    g = Phi' * r;
    [val,K] = max( abs(g) ) ;
    x(K,1) = x(K,1) + g(K,1);
    r = r - g(K,1) * Phi(:,K);
end
disp('relative error=');
disp( norm(x-x0)/norm(x) );
figure,  plot(x0,'b'); legend('原始');
figure,  plot(x,'r');  legend('重构');
figure,  plot(r);      legend('残差');