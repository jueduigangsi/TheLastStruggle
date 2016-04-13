N=256;
K=23;
M=80;
x = zeros(N,1);
q = randperm(N);
x(q(1:K)) =randn(K,1);    %ԭʼ�ź�
%% 2. �������� ���۲�ֵ���
Phi=randn(M,N); %�������� %  ��֪���󣨸�˹�ֲ���������M*N
matrixNorm = Phi.'*Phi;
matrixNorm = sqrt(diag(matrixNorm)).';
Phi = Phi./repmat(matrixNorm, [M,1]);  %ע�⣬�۲������Ҫ��һ���ģ���Ϊԭ�ӷ���Ҫ��1��
y=Phi*x ;       %������Բ���
%% 3.��MP�㷨�ع��ź�
iterations=K;                                      %  �㷨��������(m>=K)
%signal_reconstruct=zeros(1,1);                    %  ���ƽ����(��ʼֵΪ�վ���)
r_n=y;     %  �в�ֵM*1
x_rec=zeros(N,1);
for times=1:iterations
    for col=1:N                             %��֪���������������
        innerpro(col)=Phi(:,col)'*r_n;      %���������͸�֪����ÿһ�е��ڻ�
    end
    [val,pos]=max(abs(innerpro) );           %�ҳ��ڻ��о���ֵ����Ԫ�غ����Ķ�Ӧ�ĸ�֪�������pos  
    x_rec(pos)=x_rec(pos)+innerpro(pos);    %�����µĽ���x_rec
    r_n=r_n-innerpro(pos)*Phi(:,pos);       %���²в�       
end
relaerror=norm(x_rec-x)/norm(x)                           %  ����ع����
error=x_rec-x;                       %���
subplot(3,1,1);plot(x);title('origin');
subplot(3,1,2);plot(x_rec);title('reconstruct');
subplot(3,1,3);plot(error);title('�ع����');
