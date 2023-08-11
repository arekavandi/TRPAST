function [W,Z,sigma2,mu,weight]=rpast(x,data,W,Z,r,beta,sigma2,mu,Nw)
weight=0;
lambdasig=0.98;
lambdamu=0.98;

    y=W'*x;
    h=Z*y;
    g=(h/(beta+y'*h));
    e=x-W*y;
    for i=1:Nw
        A(i)=abs(norm(data(:,i)-W*W'*data(:,i),'fro')-mu)^2;
        B(i)=norm(data(:,i)-W*W'*data(:,i),'fro');
    end
    sigma2=lambdasig*sigma2+1.483*(1+5/(Nw-1))*(1-lambdasig)*median(A);
    mu=lambdamu*mu+(1-lambdamu)*median(B);
    gamma=1.96*sqrt(sigma2);
    if abs(norm(e,'fro')-mu)<gamma
    weight=1; 
    Z=(1/beta)*triu(Z-g*h')+(1/beta)*triu(Z-g*h',1)';
    W=W+e*g';
    end
end