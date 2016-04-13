clc
clear all;
close all;
M = 200;
N = 1000; 
K = 30; 
x0 = zeros(N,1);
p = randperm(N);
x0(p(1:K),1) = rand(K,1)-0.5;

Phi = randn(M,N);

% Phia = sqrt(1/M) * Phi;
% for i = 1:N
%    Phia(:,i) = Phia(:,i) / norm(Phia(:,i));
% end

y  = zeros(M,1);
y = Phi * x0;
ye2= mean(y.^2);
SNR = 40;
sgmav = sqrt( ye2*10^(-SNR/10) );
noisev = sgmav*randn(M,1);
y = y + noisev;
A = [Phi,-Phi];
b = y;
for i =1 : length(x0)
    u(i,1) = max(0,x0(i));
    v(i,1) = max(0,-x0(i));
end
x1= [ u ; v];
c= [ ones(1,length(u)) ones(1,length(v)) ];
x1= linprog(c',[],[],A,b,zeros(size(x1)),[]);
rec_x0 = x1(1:length(u)) - x1(length(u)+1:length(u)+length(v));
disp('abs_erro=');
disp(norm(rec_x0 - x0) );
disp('relative_erro=');
disp( norm(rec_x0 - x0) / norm(x0) );
stem(x0,'b*');
hold on
stem(rec_x0,'ro');
legend('原始','重构');
% figure
% stem(rec_x0 - x0,'r*');