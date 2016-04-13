function X=cosamp(Y,Compressed_Mat,m)


n=length(Y);                           
s=floor(n/4);                                         
r_n=Y;                              

sig_pos_lt=[];                        

for times=1:s                        
    
    product=abs(Compressed_Mat'*r_n);
    [val,pos]=sort(product,'descend');
    sig_pos_cr=pos(1:2*s);          
    
    sig_pos=union(sig_pos_cr,sig_pos_lt);
    
    Aug_t=Compressed_Mat(:,sig_pos);         
    
    aug_x_cr=zeros(m,1);               
    aug_x_cr(sig_pos)=(Aug_t'*Aug_t)^(-1)*Aug_t'*Y; 
    
    [val,pos]=sort(abs(aug_x_cr),'descend');
    
    X=zeros(1,m);
    X(pos(1:s))=aug_x_cr(pos(1:s));
    
    sig_pos_lt=pos(1:s);            
    
    r_n=Y-Compressed_Mat*X';
end
             



