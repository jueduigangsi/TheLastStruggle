function SP_ONE()
    clear;clc;  

    %信号长度  
    N = 256;  
    %信号满足K-稀疏  
    K = 50;  
    %观测向量的长度  
    M = floor(N*0.8);

    [y,Phi,x] = dataGen(M,N,K);  


    xp = SP_paper(Phi,y,K);  
    % [xp,support,iterCount]=SP_res(Phi,y,K);  

    figure  
    subplot(3,1,1),stem(x,'.'),xlabel('(a)原始'),axis([0,N,-1,1])  
    subplot(3,1,2),stem(xp,'.'),xlabel('(b)重构'),axis([0,N,-1,1])  
    subplot(3,1,3),stem(abs(xp-x),'.'),xlabel('(c)残差'),axis([0,N,0,max(abs(xp-x))])  

function [y,Phi,x]=dataGen(M,N,K)  
    % 产生贪婪算法所需要的数据  


    %生成-1/+1的原始信号x  
    x = zeros(N,1);  
    q = randperm(N); %y=randperm(n),是把1到n这些数随机打乱得到的一个数字序列。  
    x(q(1:K)) = sign(randn(K,1));  

    %生成测量矩阵Phi  
    Phi = randn(M,N);  
    Phi = orth(Phi')';  %求矩阵正交基,B = orth(A),B的列向量是正交向量  

    %生成观测向量  
    y = Phi*x;  
function x=SP_paper(Phi,y,K)
    %SP算法

    %获取Phi矩阵的行数和列数
    [M,N]=size(Phi);

    %初始化步骤
    %将Phi的每列与y做相关，得到一个N*1的矩阵(列向量)
    correlation=Phi'*y;
    %对correlation取绝对值后排序，按从大到小的顺序
    [var,pos] = sort(abs(correlation),'descend');
    %声明一个空集T，用于记录Phi的列数标值
    T=[];T=union(T,pos(1:K));
    y_r=resid_paper(y,Phi(:,T));

    %迭代
    %使用如下形式的do---while结构
    % while(1)
    %     if(condition)
    %         break;
    %     end
    % end
    count=1;
    while(1)
        %根据残差计算待增加的列数，得到T_add
        correlation=Phi'*y_r;
        [var,pos] = sort(abs(correlation),'descend');
        T_add=union([],pos(1:K));

        %合并已有的T和T_add
        T=union(T,T_add);

        %
        x_p=((Phi(:,T)'*Phi(:,T))\eye(length(T)))*Phi(:,T)'*y;%proj_paper(y,Phi(:,T));
        %更新下标记录T
        [var,pos] = sort(abs(x_p),'descend');
        %取前K个最大值
        T=union([],T(pos(1:K)));

        %计算新的残差
        y_r_n=resid_paper(y,Phi(:,T));

        %判断是否退出循环,且置为最大迭代100次
        if(norm(y_r_n)>=norm(y_r) || count>100)
            break;
        end

        %若不退出循环,进行新一轮的迭代
        y_r=y_r_n;

        count=count+1;
    end

    %退出循环后，做最后的数据输出
    x=zeros(N,1);
    x(T)=((Phi(:,T)'*Phi(:,T))\eye(length(T)))*Phi(:,T)'*y;








function y_r=resid_paper(y,Phi)
    %计算y在Phi上的投影残差

    %获取矩阵Phi的行数和列数,M没有用
    [M,N]=size(Phi);

    %判断矩阵(Phi'*Phi)是否可逆
    if(rank(Phi'*Phi)~=N)
        error('矩阵不可逆');
    end

    y_p=Phi*((Phi'*Phi)\eye(N))*Phi'*y;
    y_r=y-y_p;

