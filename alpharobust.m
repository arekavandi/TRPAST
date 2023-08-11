function [W,Z,weight]=alpharobust(x,W,Z,r,beta,alpha);

    y=W'*x; 
    delta=sum((x-W*y).^2);
    weight=exp(((alpha-1)/2)*delta)+0.000001;
       
    h=Z*y;
    g=(h*weight/(beta+y'*h*weight));
    Z=1/beta*(Z-g*h');
    e=x-W*y;
    
    taw=norm(g)^(-2)*((1+norm(e)^2*norm(g)^2)^(-0.5)-1);
    ep=taw*W*g+(1+taw*norm(g)^2)*e;
    W=W+ep*g';
    
end