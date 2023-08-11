function [W,Z]=past(x,W,Z,r,beta)
    y=W'*x;
    h=Z*y;
    g=(h/(beta+y'*h));
    Z=(1/beta)*triu(Z-g*h')+(1/beta)*triu(Z-g*h',1)';
    e=(x-W*y);
    W=W+e*g';
end