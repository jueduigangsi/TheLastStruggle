
clc;clear

%%  1. 时域测试信号生成
% K=7;      %  稀疏度(做FFT可以看出来)
% N=256;    %  信号长度
% M=64;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
% f1=50;    %  信号频率1
% f2=100;   %  信号频率2
% f3=200;   %  信号频率3
% f4=400;   %  信号频率4
% fs=800;   %  采样频率
% ts=1/fs;  %  采样间隔
% Ts=1:N;   %  采样序列
% x=cos(2*pi*Ts*ts);  %  完整信号

M=50;
K=20;
k = 1:100;
N = 100;
x = sin(4*k*pi/N);
%%  2.  时域信号压缩传感
Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
s=Phi*x.';                                        %  获得线性测量 

%%  3.  正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
m=2*K;                                            %  算法迭代次数(m>=K)
Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)

hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=s;                                            %  残差值

for times=1:m;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小  aug_y就是估计的信号！！
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
hat_y(pos_array)=aug_y;                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.');                         %  做逆傅里叶变换重构得到时域信号

%%  4.  恢复信号和原始信号对比
figure(1);
hold on;
subplot(3,1,1);
 stem(x,'.');                           
                     %  原始信号
 xlabel('(a) 原始');
 subplot(3,1,2);
  stem(hat_x,'.');                  %  重建信号

 xlabel('(b) 重构');
 subplot(3,1,3);
 plot(hat_x-x');          
 xlabel('(c) 误差');
norm(hat_x.'-x)/norm(x)                           %  重构误差