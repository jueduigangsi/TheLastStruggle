clc
clear all;
close all;
M = 800;   %�۲��źų���
N = 256*256;  %ϡ���źų���
K = 23;   %ϡ���
%% -----1.����ϡ���ΪK��ϡ���ź�-----
x0 = zeros(N,1);
p = randperm(N);
x0(p(1:K),1) = randn(K,1);
%------  ��˹��֪����Phi   -------------
Phi = sqrt(1/M) * randn(M,N);
for i = 1:N
    Phi(:,i) = Phi(:,i) / norm(Phi(:,i));
end
%-------- �������� y  ----------
y  = zeros(M,1);
y = Phi * x0;
%% -----2. MP Reconstruction ------------
x = zeros(N,1);       %x0�ıƽ��ź�x
times =  K;           %�������� = ϡ���
g = zeros(N,1);       %�����͸�֪�����ڻ�
r  = y;               %������ʼ��Ϊy
for n=0:times
    g = Phi' * r;
    [val,K] = max( abs(g) ) ;
    x(K,1) = x(K,1) + g(K,1);
    r = r - g(K,1) * Phi(:,K);
end
disp('relative error=');
disp( norm(x-x0)/norm(x) );
figure,  plot(x0,'b'); legend('ԭʼ');
figure,  plot(x,'r');  legend('�ع�');
figure,  plot(r);      legend('�в�');