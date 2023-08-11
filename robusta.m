function [W,Z,K,weight]=robusta(x,K,W,Z,r,beta);
    u=K*x;
    delta=x'*u;
    v=u/(beta+delta);
    K=1/beta*(K-u*v');
    weight=1/delta;
    
    y=W'*x;
    h=Z*y;
    g=(h*weight/(beta+y'*h*weight));
    Z=1/beta*(Z-g*h');
    e=x-W*y;
    
    taw=norm(g)^(-2)*((1+norm(e)^2*norm(g)^2)^(-0.5)-1);
    ep=taw*W*g+(1+taw*norm(g)^2)*e;
    W=W+ep*g';
    
end