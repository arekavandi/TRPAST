function [W,Z]=opast(x,W,Z,r,beta)
    y=W'*x;
    q=(1/beta)*Z*y;
    gamma=(1/(1+y'*q));
    Z=(1/beta)*Z-gamma*q*q';
    p=gamma*(x-W*y);
    taw=(1/(norm(q)^2))*(1/sqrt(1+(norm(p)^2)*(norm(q)^2))-1);
    pprime=taw*W*q+(1+taw*norm(q)^2)*p;
    W=W+pprime*q';
end