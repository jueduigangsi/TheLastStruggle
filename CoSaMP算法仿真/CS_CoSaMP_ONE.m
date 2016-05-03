 clear;clc;  

    %�źų���  
    N = 256;  
    %�ź�����K-ϡ��  
    K = 50;  
    %�۲������ĳ���  
    M = floor(N*0.8);
    
    x = zeros(N,1);  
    q = randperm(N); %y=randperm(n),�ǰ�1��n��Щ��������ҵõ���һ���������С�  
    x(q(1:K)) = sign(randn(K,1));  

    %���ɲ�������Phi  
    Phi = randn(M,N);  
    Phi = orth(Phi')';  %�����������,B = orth(A),B������������������  

    %���ɹ۲�����  
    y = Phi*x; 
    xp=cosamp(y,Phi,N);
    xp=xp';
   figure  
    subplot(3,1,1),stem(x,'.'),xlabel('(a)ԭʼ'),axis([0,N,-1,1])  
    subplot(3,1,2),stem(xp,'.'),xlabel('(b)�ع�'),axis([0,N,-1,1])  
    subplot(3,1,3),stem(abs(xp-x),'.'),xlabel('(c)�в�'),axis([0,N,0,max(abs(xp-x))])  