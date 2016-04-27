%一维信号BP重构算法
clc
clear all;
close all;
M = 500;   %观测信号长度
N = 1000;  %稀疏信号长度
K = 50;   %稀疏度
%% -----1.生成稀疏度为K的稀疏信号-----
x0 = zeros(N,1);
p = randperm(N);
x0(p(1:K),1) = rand(K,1)-0.5;
%------  高斯感知矩阵Phi   -------------
Phi = sqrt(1/M) * randn(M,N);
for i = 1:N
    Phi(:,i) = Phi(:,i) / norm(Phi(:,i));
end
%-------- 测量向量 y  ----------
y  = zeros(M,1);
y = Phi * x0;
%%  -----2. 含高斯白噪声观测 SNR=40dB --------
ye2 = mean(y.^2);
SNR = 40;
sgmav = sqrt( ye2*10^(-SNR/10) );
noisev = sgmav*randn(M,1);
y = y + noisev;
%% ------3.BP重构算法 ----------------
A = [Phi,-Phi];
b = y;
for i =1 : length(x0)
    u(i,1) = max(0,x0(i));
    v(i,1) = max(0,-x0(i));
end
x = [ u ; v];
c = [ ones(1,length(u)) ones(1,length(v)) ];
x = linprog(c',[],[],A,b,zeros(size(x)),[]);
rec_x0 = x(1:length(u)) - x(length(u)+1:length(u)+length(v));
%% ----4.重建质量衡量--------------
disp('abs_err=');
disp(norm(rec_x0 - x0) );
disp('relative_erro=');
disp( norm(rec_x0 - x0) / norm(x0) );
SNR = 20 * log10( norm(x0)/norm(rec_x0-x0) )
stem(x0,'.');
figure;
stem(rec_x0,'.');
legend('原始','重构');
figure
 stem(rec_x0 - x0,'.');