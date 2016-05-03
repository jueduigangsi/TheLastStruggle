 clear;clc;  

    %信号长度  
    N = 256;  
    %信号满足K-稀疏  
    K = 50;  
    %观测向量的长度  
    M = floor(N*0.8);
    
    x = zeros(N,1);  
    q = randperm(N); %y=randperm(n),是把1到n这些数随机打乱得到的一个数字序列。  
    x(q(1:K)) = sign(randn(K,1));  

    %生成测量矩阵Phi  
    Phi = randn(M,N);  
    Phi = orth(Phi')';  %求矩阵正交基,B = orth(A),B的列向量是正交向量  

    %生成观测向量  
    y = Phi*x; 
    xp=cosamp(y,Phi,N);
    xp=xp';
   figure  
    subplot(3,1,1),stem(x,'.'),xlabel('(a)原始'),axis([0,N,-1,1])  
    subplot(3,1,2),stem(xp,'.'),xlabel('(b)重构'),axis([0,N,-1,1])  
    subplot(3,1,3),stem(abs(xp-x),'.'),xlabel('(c)残差'),axis([0,N,0,max(abs(xp-x))])  