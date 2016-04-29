clc
close all
clear all
I = imread('lena128.bmp');
[m,n] = size(I);
N = m * n;  %ϡ���źų���
vectorI = reshape(I,N,1);
x0 = dct2(vectorI);
sampleRate = 0.3;  %������ = M/N
M = round( N * sampleRate );
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
times =  M/2        %�������� = ϡ���
g = zeros(N,1);       %�����͸�֪�����ڻ�
r  = y;               %������ʼ��Ϊy
for count=0:times
    count;
    g = Phi' * r;
    [val,K] = max( abs(g) ) ;
    x(K,1) = x(K,1) + g(K,1);
    r = r - g(K,1) * Phi(:,K);
end
idct_x = idct(x);
rec_x = reshape(idct_x,m,n);
%% --------3.�ؽ���������--------------------
relative_err = norm(x-x0) / norm(x)
MSE = norm(x-x0) / (m*n)
SNR = 10*log10( norm(x)^2 / MSE )
Mat_rate = 1 - norm(abs(x)-abs(x0))/norm(abs(x)+abs(x0))
figure
subplot(121),imshow(I,[]),title('orignal image');
subplot(122),imshow(rec_x,[]), title('reconstruct image');


%%  ----�������ݼ�¼--------------------------
% times = M,sampleRate = 0.5;
% relative_err =
%     0.1278
% MSE =
%     0.2658
% SNR =
%    84.3585

% times = M/2,sampleRate = 0.3;
% ��������Ϊ614���ؽ��ź�ϡ���K=588
% relative_err =
%     0.1593
% MSE =
%     0.3299
% SNR =
%    83.3882
% Mat_rate =
%     0.9235
