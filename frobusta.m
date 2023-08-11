function [W,Z,Tx,Ty,weight]=frobusta(x,Tx,Ty,W,Z,r,beta); 
   [dim d2]=size(x);
    y=W'*x;
    Tx=beta*Tx+norm(x)^2;
    Ty=beta*Ty+norm(y)^2;
    h=Z*y;
    delta=y'*h+((norm(x)^2-norm(y)^2)*(dim-r))/(Tx-Ty);
    weight=1/delta;
    
    g=(h*weight/(beta+y'*h*weight));
    Z=(1/beta)*(Z-g*h');
    e=x-W*y;
    
    taw=norm(g)^(-2)*((1+norm(e)^2*norm(g)^2)^(-0.5)-1);
    ep=taw*W*g+(1+norm(g)^2)*e;
    W=W+ep*g';
    
end