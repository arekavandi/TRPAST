function [W, Z,K,Tx,Ty,w] = new_RFAPI(X,W,Z,K,beta,Tx,Ty, r)

% X is the n x N data matrix collecting the N observation (n x 1) vectors 
% beta is the forgetting factor
% alpha is the tuning factor 
% r is the number of principal eigenvectors
% V (n x r):  the exact principal subspace basis used for the subspace error evaluation

%New Robust Fapi Processing

    y = W'*X;
    h = Z*y;
    u=K*X;
    eps = (norm(X))^2 - (norm(y))^2;
%     delta=y'*h+(eps)*(size(W,1)-size(W,2))/(Tx-Ty+0.001);
   
  
    delta=X'*u;
    w=1/delta;
    v=u/(beta+delta);
    K=1/beta*(K-u*v');
%     Tx=beta*Tx+w*norm(X)^2;
%     Ty=beta*Ty+w*norm(y)^2;
    g = h*w/(beta + y'*h*w);
    gn = (norm(g))^2;
    %FAPI main section
    tau = eps/(1+ eps*gn + sqrt(1 + eps*gn) );
    eta = 1 - tau*gn;
    yp = eta*y + tau*g;
    hp = Z'* yp;
    e = (tau/eta)*(Z*g - (hp'*g)*g);
    Z = (1/beta)* (Z - g*hp' + e*g');
    ep = eta* X - W*yp;
    W = W + ep*g';

end

    
