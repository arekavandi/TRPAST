clc
clear all
close all

montecarlorun=50;

for index=1:5
    if index==5
        index=20;
    end
        if index==2
        index=1.02;
    end
for it=1:montecarlorun
N=5000;
dim=20;
r=2;
beta=0.999;
ASNR=10;
alpha=0.999-(index-1)*0.2;
sigma=0.2;
ep=0.2;


A=rand(dim,dim)-0.5;
% B=A'*A;
 [U S V]=svd(A);
H=A(:,1:r);

s=H*rand(r,100);
norms=mean(sum(s.^2));
scale=sqrt((sigma^2*dim*10^(ASNR/10))/norms);

% for i=1:r
%     H(:,i)=H(:,i)/norm(H(:,i));
% end

for i=1:N
    a=rand(r,1);
    theta=scale*a;
%     X(:,i)=H*theta+wblrnd(0.5,0.5,[dim 1]);
     X(:,i)=H*theta+random('normal',0,sigma,dim,1); 
end


 %outliers vector level
 outlier=ones(dim,N);
 for i=1:N
         if rand<ep
             outlier(:,i)=zeros(dim,1);
            X(:,i)= X(:,i)+random('normal',0,10,dim,1);
         end
 end


% 
% % Outliers element level
% outlier=ones(dim,N);
% for i=1:N
%     for j=1:dim
%         if rand<ep
%             outlier(j,i)=0;
%            X(j,i)=X(j,i)+random('normal',5,0.1);
%         end
%     end
% end


%%% Subspace Estimation %%%%
C=eye(dim);
K=eye(dim);
W=[eye(r);rand(dim-r,r)*0];
Z=eye(r);

Wapi=W;
Zapi=Z;

Wfapi=W;
Zfapi=Z;

Wrfapi=W;
Zrfapi=Z;
Tx=0;
Ty=0;
Krfapi=K;

Wrobusta=W;
Zrobusta=Z;
Krobusta=K;

Walpha=W;
Zalpha=Z;

Walpha2=W;
for i=1:dim
    Zalpha2(:,:,i)=Z;
end


Txfrobusta=dim;
Tyfrobusta=r;
Wfrobusta=W;
Zfrobusta=Z;

Wpast=W;
Zpast=Z;

Wopast=W;
Zopast=Z;

murpast=0;
sigma2rpast=1;
Wrpast=W;
Zrpast=Z;

B=H*(H'*H)^(-0.5);

weightalpha2=ones(dim,1);
weightrpast=zeros(dim,N);
% mu=Walpha2'*X(:,1);
for i=1:N
%     [Wapi,Zapi]=api(X(:,i),Wapi,Zapi,r,beta);
%     lossapi(i,it)=trace(Wapi'*(eye(dim)-B*B')*Wapi)/trace(Wapi'*(B*B')*Wapi);
%     OEapi(i,it)=norm(Wapi'*Wapi-eye(r),'fro');
    
%     [Wfapi,Zfapi]=fapi(X(:,i),Wfapi,Zfapi,r,beta);
%     lossfapi(i,it)=trace(Wfapi'*(eye(dim)-B*B')*Wfapi)/trace(Wfapi'*(B*B')*Wfapi);
%     OEfapi(i,it)=norm(Wfapi'*Wfapi-eye(r),'fro');
    
%     [Wrfapi,Zrfapi,Krfapi, Tx,Ty,wei]=new_RFAPI(X(:,i),Wrfapi,Zrfapi,Krfapi,beta, Tx,Ty, r);
%     weirfapi(:,i)=wei*ones(dim,1);
%     lossrfapi(i,it)=trace(Wrfapi'*(eye(dim)-B*B')*Wrfapi)/trace(Wrfapi'*(B*B')*Wrfapi);
%     OErfapi(i,it)=norm(Wrfapi'*Wrfapi-eye(r),'fro');
    
    [Wrobusta,Zrobusta,Krobusta,e]=robusta(X(:,i),Krobusta,Wrobusta,Zrobusta,r,beta);
    weight(:,i)=e*ones(dim,1);
    lossrobusta(i,it)=trace(Wrobusta'*(eye(dim)-B*B')*Wrobusta)/trace(Wrobusta'*(B*B')*Wrobusta);
    OErobusta(i,it)=norm(Wrobusta'*Wrobusta-eye(r),'fro');
    
%     [Wfrobusta,Zfrobusta,Txfrobusta,Tyfrobusta,weight1(i)]=frobusta(X(:,i),Txfrobusta,Tyfrobusta,Wfrobusta,Zfrobusta,r,beta);
%     lossfrobusta(i,it)=trace(Wfrobusta'*(eye(dim)-B*B')*Wfrobusta)/trace(Wfrobusta'*(B*B')*Wfrobusta);
%     OEfrobusta(i,it)=norm(Wfrobusta'*Wfrobusta-eye(r),'fro');
    
%     [Wpast,Zpast]=past(X(:,i),Wpast,Zpast,r,beta);
%     losspast(i,it)=trace(Wpast'*(eye(dim)-B*B')*Wpast)/trace(Wpast'*(B*B')*Wpast);
%     OEpast(i,it)=norm(Wpast'*Wpast-eye(r),'fro');
    
%     [Wopast,Zopast]=opast(X(:,i),Wopast,Zopast,r,beta);
%     lossopast(i,it)=trace(Wopast'*(eye(dim)-B*B')*Wopast)/trace(Wopast'*(B*B')*Wopast);
%     OEopast(i,it)=norm(Wopast'*Wopast-eye(r),'fro');
    
    [Walpha,Zalpha,e]=alpharobust(X(:,i),Walpha,Zalpha,r,beta,alpha);
    weightalpha(:,i)=e*ones(dim,1);
    lossalpharobust(i,it)=trace(Walpha'*(eye(dim)-B*B')*Walpha)/trace(Walpha'*(B*B')*Walpha);
    OEalpharobust(i,it)=norm(Walpha'*Walpha-eye(r),'fro');
    
    
% % weightalpha2
% 
%     [Walpha2,Zalpha2,e,alpha2]=alpharobust2(X(:,i),Walpha2,Zalpha2,r,beta,alpha2,weightalpha2,ep,i);
%     weightalpha2(:,i)=e;
% %     Walpha2=Walpha2*(Walpha2'*Walpha2)^(-0.5);
%     lossalpharobust2(i,it)=trace(Walpha2'*(eye(dim)-B*B')*Walpha2)/trace(Walpha2'*(B*B')*Walpha2);
%     OEalpharobust2(i,it)=norm(Walpha2'*Walpha2-eye(r),'fro'); 

% % 
% % alpha2
%     if i>=Nw
%     [Wrpast,Zrpast,sigma2rpast,murpast,e]=rpast(X(:,i),X(:,i-Nw+1:i),Wrpast,Zrpast,r,beta,sigma2rpast,murpast,Nw);
%     weightrpast(:,i)=e*ones(dim,1);
%     lossrpast(i,it)=trace(Wrpast'*(eye(dim)-B*B')*Wrpast)/trace(Wrpast'*(B*B')*Wrpast);
%     OErpast(i,it)=norm(Wrpast'*Wrpast-eye(r),'fro');
%     end
%     
% 
% if i==20
%     gg
% end
end

end

 
    

if index==1
    figure
    semilogy(mean(lossrobusta'),'--','LineWidth',2)
    hold on
    semilogy(mean(lossalpharobust'),'LineWidth',2)
    hold on

else
    semilogy(mean(lossalpharobust'),'LineWidth',2)
    hold on
end
end
    xlabel('Time', 'Interpreter', 'LaTeX')
    ylabel('SEP', 'Interpreter', 'LaTeX')
    grid on
    legend({'ROBUSTA@2018','A1, alpha=1','A1, alpha=0.99','A1, alpha=0.6','A1, alpha=0.4','A1, alpha=-3'}, ...
        'Interpreter', 'LaTeX')
%print('Performance-H1085','-depsc')
