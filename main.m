
%This code was developed by Aref Miri Rekavandi for the paper: Rekavandi, A. M., Seghouane, A. K., &
%Abed-Meraim, K. (2023). TRPAST: A tunable and robust projection approximation subspace tracking method.
%IEEE Transactions on Signal Processing. 
%If you use this code in your study, kindly cite the aforementioned paper.

clc
clear all
close all
%% Parameters %%
montecarlorun=5;       %Number of trials
N=10000;                %Number of samples in each trial
dim=20;                 %Signal dimension
r=2;                    %Subspace rank
beta=0.999;             %Forgetting factor
ASNR=10;                %SNR
alpha=0.5;              %Alpha in alpha-divergence
alpha2=0.5*ones(dim,1);
sigma=0.2;              %Noise std
ep=0.2;                 %Contamination portion
Nw=10;                  %Window size RPAST algorithm
%%
for it=1:montecarlorun
fprintf('Trial %d out of %d\n',it,montecarlorun)
A=rand(dim,dim)-0.5;
% B=A'*A;
 [U S V]=svd(A);
H=A(:,1:r);

s=H*rand(r,100);
norms=mean(sum(s.^2));
scale=sqrt((sigma^2*dim*10^(ASNR/10))/norms);


for i=1:N
    a=rand(r,1);
    theta=scale*a;
     X(:,i)=H*theta+random('normal',0,sigma,dim,1); 
end


%  %outliers vector level
%  outlier=ones(dim,N);
%  for i=1:N
%          if rand<ep
%              outlier(:,i)=zeros(dim,1);
%             X(:,i)= X(:,i)+random('normal',0,10,dim,1);
%          end
%  end



% Outliers element level
outlier=ones(dim,N);
for i=1:N
    for j=1:dim
        if rand<ep
            outlier(j,i)=0;
           X(j,i)=X(j,i)+random('normal',5,0.1);
        end
    end
end


%%% Subspace Estimation %%%%
C=eye(dim);
K=eye(dim);
W=[eye(r);rand(dim-r,r)*0];
Z=eye(r);

Wapi=W; Zapi=Z;   

Wfapi=W;Zfapi=Z;    

Wrfapi=W; Zrfapi=Z; Tx=0; Ty=0; Krfapi=K;

Wrobusta=W;Zrobusta=Z;Krobusta=K;     

Walpha=W;Zalpha=Z;    Walpha2=W;

for i=1:dim
    Zalpha2(:,:,i)=Z;
end


Txfrobusta=dim;  Tyfrobusta=r;  Wfrobusta=W; Zfrobusta=Z;

Wpast=W; Zpast=Z;

Wopast=W; Zopast=Z;

murpast=0; sigma2rpast=1; Wrpast=W; Zrpast=Z; 

B=H*(H'*H)^(-0.5);

weightalpha2=ones(dim,1);
weightrpast=zeros(dim,N);

for i=1:N
    [Wapi,Zapi]=api(X(:,i),Wapi,Zapi,r,beta);
    lossapi(i,it)=trace(Wapi'*(eye(dim)-B*B')*Wapi)/trace(Wapi'*(B*B')*Wapi);
    OEapi(i,it)=norm(Wapi'*Wapi-eye(r),'fro');
    
    [Wfapi,Zfapi]=fapi(X(:,i),Wfapi,Zfapi,r,beta);
    lossfapi(i,it)=trace(Wfapi'*(eye(dim)-B*B')*Wfapi)/trace(Wfapi'*(B*B')*Wfapi);
    OEfapi(i,it)=norm(Wfapi'*Wfapi-eye(r),'fro');
    
    [Wrfapi,Zrfapi,Krfapi, Tx,Ty,wei]=new_RFAPI(X(:,i),Wrfapi,Zrfapi,Krfapi,beta, Tx,Ty, r);
    weirfapi(:,i)=wei*ones(dim,1);
    lossrfapi(i,it)=trace(Wrfapi'*(eye(dim)-B*B')*Wrfapi)/trace(Wrfapi'*(B*B')*Wrfapi);
    OErfapi(i,it)=norm(Wrfapi'*Wrfapi-eye(r),'fro');
    
    [Wrobusta,Zrobusta,Krobusta,e]=robusta(X(:,i),Krobusta,Wrobusta,Zrobusta,r,beta);
    weight(:,i)=e*ones(dim,1);
    lossrobusta(i,it)=trace(Wrobusta'*(eye(dim)-B*B')*Wrobusta)/trace(Wrobusta'*(B*B')*Wrobusta);
    OErobusta(i,it)=norm(Wrobusta'*Wrobusta-eye(r),'fro');
    
    [Wfrobusta,Zfrobusta,Txfrobusta,Tyfrobusta,weight1(i)]=frobusta(X(:,i),Txfrobusta,Tyfrobusta,Wfrobusta,Zfrobusta,r,beta);
    lossfrobusta(i,it)=trace(Wfrobusta'*(eye(dim)-B*B')*Wfrobusta)/trace(Wfrobusta'*(B*B')*Wfrobusta);
    OEfrobusta(i,it)=norm(Wfrobusta'*Wfrobusta-eye(r),'fro');
    
    [Wpast,Zpast]=past(X(:,i),Wpast,Zpast,r,beta);
    losspast(i,it)=trace(Wpast'*(eye(dim)-B*B')*Wpast)/trace(Wpast'*(B*B')*Wpast);
    OEpast(i,it)=norm(Wpast'*Wpast-eye(r),'fro');
    
    [Wopast,Zopast]=opast(X(:,i),Wopast,Zopast,r,beta);
    lossopast(i,it)=trace(Wopast'*(eye(dim)-B*B')*Wopast)/trace(Wopast'*(B*B')*Wopast);
    OEopast(i,it)=norm(Wopast'*Wopast-eye(r),'fro');
    
    %First proposed method
    [Walpha,Zalpha,e]=alpharobust(X(:,i),Walpha,Zalpha,r,beta,alpha);
    weightalpha(:,i)=e*ones(dim,1);
    lossalpharobust(i,it)=trace(Walpha'*(eye(dim)-B*B')*Walpha)/trace(Walpha'*(B*B')*Walpha);
    OEalpharobust(i,it)=norm(Walpha'*Walpha-eye(r),'fro');
    
    
    %Second Proposed method
    [Walpha2,Zalpha2,e,alpha2]=alpharobust2(X(:,i),Walpha2,Zalpha2,r,beta,alpha2,weightalpha2,ep,i);
    weightalpha2(:,i)=e;
    lossalpharobust2(i,it)=trace(Walpha2'*(eye(dim)-B*B')*Walpha2)/trace(Walpha2'*(B*B')*Walpha2);
    OEalpharobust2(i,it)=norm(Walpha2'*Walpha2-eye(r),'fro'); 

    if i>=Nw
    [Wrpast,Zrpast,sigma2rpast,murpast,e]=rpast(X(:,i),X(:,i-Nw+1:i),Wrpast,Zrpast,r,beta,sigma2rpast,murpast,Nw);
    weightrpast(:,i)=e*ones(dim,1);
    lossrpast(i,it)=trace(Wrpast'*(eye(dim)-B*B')*Wrpast)/trace(Wrpast'*(B*B')*Wrpast);
    OErpast(i,it)=norm(Wrpast'*Wrpast-eye(r),'fro');
    end
    
end

end
      
     
start=9850;
last=10000;
start1=1;
last1=150;
figure
subplot(6,1,1)
imshow([outlier(:,start1:last1) outlier(:,start:last)])
title({'Ideal Weight'}, ...
        'Interpreter', 'LaTeX')
    colorbar
  subplot(6,1,5)
 imshow(log10([weight(:,start1:last1) weight(:,start:last)])/max(max(log10(weight(:,start:last)))),[])
 title({'Weight ROBUSTA@2018'}, ...
        'Interpreter', 'LaTeX') 
    colorbar
      subplot(6,1,4)
 imshow([weightrpast(:,start1:last1) weightrpast(:,start:last)],[])
 title({'Weight RPAST@2005'}, ...
        'Interpreter', 'LaTeX') 
    colorbar
 subplot(6,1,2)
 imshow([weightalpha(:,start1:last1) weightalpha(:,start:last)])
 title({'Weight A1'}, ...
        'Interpreter', 'LaTeX')
    colorbar
 subplot(6,1,3)
 imshow([weightalpha2(:,start1:last1) weightalpha2(:,start:last)])
 title({'Weight A2'}, ...
        'Interpreter', 'LaTeX')
    colorbar
 subplot(6,1,6)
 imshow(log10([weirfapi(:,start1:last1) weirfapi(:,start:last)])/max(max(log10(weirfapi(:,start:last)))),[])
 title({'Weight MFAPI@2022'}, ...
        'Interpreter', 'LaTeX')
    colorbar
    
    
figure
subplot(1,2,1)

semilogy(mean(losspast'),'LineWidth',2)
hold on
semilogy(mean(lossopast'),'LineWidth',2)
hold on
semilogy(mean(lossfapi'),'LineWidth',2)
hold on
semilogy(mean(lossrfapi'),'LineWidth',2)
hold on
semilogy(mean(lossrpast'),'LineWidth',2)
hold on
semilogy(mean(lossrobusta'),'--','LineWidth',2)

hold on
semilogy(mean(lossalpharobust'),'LineWidth',2)
hold on
semilogy(mean(lossalpharobust2'),'k','LineWidth',2)
    xlabel('Time', 'Interpreter', 'LaTeX')
    ylabel('SEP', 'Interpreter', 'LaTeX')
    grid on
    legend({'PAST@1995','OPAST@2000','FAPI@2005','MFAPI@2022','RPAST@2005','ROBUSTA@2018','A1','A2'}, ...
        'Interpreter', 'LaTeX')
subplot(1,2,2)

semilogy(mean(OEpast'),'LineWidth',2)
hold on
semilogy(mean(OEopast'),'LineWidth',2)
hold on
semilogy(mean(OEfapi'),'LineWidth',2)
hold on
semilogy(mean(OErfapi'),'LineWidth',2)
hold on
semilogy(mean(OErpast'),'LineWidth',2)
hold on
semilogy(mean(OErobusta'),'LineWidth',2)

hold on
semilogy(mean(OEalpharobust'),'LineWidth',2)
hold on
semilogy(mean(OEalpharobust2'),'k','LineWidth',2)
    xlabel('Time', 'Interpreter', 'LaTeX')
    ylabel('OE', 'Interpreter', 'LaTeX')
    grid on
    legend({'PAST@1995','OPAST@2000','FAPI@2005','MFAPI@2022','RPAST@2005','ROBUSTA@2018','A1','A2'}, ...
        'Interpreter', 'LaTeX')

