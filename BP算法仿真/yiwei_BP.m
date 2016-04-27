%һά�ź�BP�ع��㷨
clc
clear all;
close all;
M = 500;   %�۲��źų���
N = 1000;  %ϡ���źų���
K = 50;   %ϡ���
%% -----1.����ϡ���ΪK��ϡ���ź�-----
x0 = zeros(N,1);
p = randperm(N);
x0(p(1:K),1) = rand(K,1)-0.5;
%------  ��˹��֪����Phi   -------------
Phi = sqrt(1/M) * randn(M,N);
for i = 1:N
    Phi(:,i) = Phi(:,i) / norm(Phi(:,i));
end
%-------- �������� y  ----------
y  = zeros(M,1);
y = Phi * x0;
%%  -----2. ����˹�������۲� SNR=40dB --------
ye2 = mean(y.^2);
SNR = 40;
sgmav = sqrt( ye2*10^(-SNR/10) );
noisev = sgmav*randn(M,1);
y = y + noisev;
%% ------3.BP�ع��㷨 ----------------
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
%% ----4.�ؽ���������--------------
disp('abs_err=');
disp(norm(rec_x0 - x0) );
disp('relative_erro=');
disp( norm(rec_x0 - x0) / norm(x0) );
SNR = 20 * log10( norm(x0)/norm(rec_x0-x0) )
stem(x0,'.');
figure;
stem(rec_x0,'.');
legend('ԭʼ','�ع�');
figure
 stem(rec_x0 - x0,'.');