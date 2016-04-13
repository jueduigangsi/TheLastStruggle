function X=omp(Y,Compressed_Mat,m)

n=length(Y);
s=floor(n/4);                                     %  ����ֵά��
X=zeros(1,m);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=Y;                                            %  �в�ֵ 

for times=1:s;                                  %  ��������(ϡ����ǲ�����1/4)

    product=abs(Compressed_Mat'*r_n);
    
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,Compressed_Mat(:,pos)];                   %  ��������
    Compressed_Mat(:,pos)=zeros(n,1);                      %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼򵥽������㣩
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*Y;           %  ��С����,ʹ�в���С
    r_n=Y-Aug_t*aug_x;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
    
end
X(pos_array)=aug_x;                           %  �ع������� 


