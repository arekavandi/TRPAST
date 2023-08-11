function [W,Z,out,alpha]=alpharobust2(x,W,Z,r,beta,alpha,weightalpha,ep,num);
    dim=size(x,1);
    
    y=W'*x; 
        
    for j=1:dim
        
        if mod(num,25)==10 && num>100
        if (sum(weightalpha(j,num-100:num-1)<0.05)/100)>(ep)
            alpha(j)=alpha(j)+0.01;
        elseif (sum(weightalpha(j,num-100:num-1)<0.05)/100)<(ep-0.01)
            alpha(j)=alpha(j)-0.01;
        end
        if (sum(weightalpha(j,num-100:num-1)>0.5)/100)<(1-ep)
            alpha(j)=alpha(j)+0.005;
        end
        if alpha(j)>1
            alpha(j)=1;
        end
        end
    deltaj=(x(j))^2;
    weight(j)=exp(((alpha(j)-1)/2)*deltaj)+0.000000001;
       
    h=Z(:,:,j)*y;
    gj=(h*weight(j)/(beta+y'*h*weight(j)));
    Z(:,:,j)=1/beta*(Z(:,:,j)-gj*h');
    ej=x(j)-W(j,:)*y;
    
    W(j,:)=W(j,:)+ej*gj';
    end
    W=W*(W'*W)^(-0.5);
    out=weight';
%     mu=beta*mu+W'*(weight'.*x);
    
end