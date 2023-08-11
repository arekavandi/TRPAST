function [W,Z]=fapi(x,W,Z,r,beta)
    y=W'*x;
    h=Z*y;
    g=(h/(beta+y'*h));
    
    epsilon2=norm(x)^2-norm(y)^2;
    taw=epsilon2/(1+epsilon2*norm(g)^2+sqrt(1+epsilon2*norm(g)^2));
    eta=1-taw*norm(g)^2;
    yp=eta*y+taw*g;
    hp=Z'*yp;
    e=(taw/eta)*(Z*g-hp'*g*g);
    Z=1/beta*(Z-g*hp'+e*g');
    ep=eta*x-W*yp;
    W=W+ep*g';
end