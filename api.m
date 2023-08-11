function [W,Z]=api(x,W,Z,r,beta)
    y=W'*x;
    h=Z*y;
    g=(h/(beta+y'*h));
    e=(x-W*y);
    Capt=(eye(r)+norm(e)^2*g*g')^(-0.5);
    Z=1/beta*Capt'*(eye(r)-g*y')*Z*(Capt')^(-1);
    W=(W+e*g')*Capt;
end