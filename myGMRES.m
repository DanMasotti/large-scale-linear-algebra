function [x,res] = myGMRES(A,b,x0,m,be)

r0 = b-A*x0;
beta = norm(r0,2);

n = size(A,1);
V = zeros(n,m);
H = zeros(m,m);

V(:,1) = r0/beta;

tol = 1e-12;

% for j=1:m
% res(1) = 100000;
j = 1;
% while (res(j) >= tol)
for j=1:m
    w = A*V(:,j);
    
    for i=1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    
    H(j+1,j) = norm(w);
    
    if (H(j+1,j) == 0)
        disp('hello');
        break;
    else
        V(:,j+1) = w/H(j+1,j);
    end
    
    
    e1j = zeros(j+1,1);
    e1j(1) = 1;
    
    Hj = H(1:j+1,1:j);
%     y = (Hj'*Hj)\(Hj(1,:)'*beta);
    
%     [Q,R] = qr(Hj);
%     y = R\(Q(1,:)'*beta);
    
    y = Hj\(beta*e1j);
    x = x0 + V(:,1:j)*y;
    
    res(j) = norm(b-A*x,2);
%     res(j) = norm(beta*e1j - Hj*y,2);%/norm(be);
    
    
%     y = H(1:j,1:j)\(beta*e1j);
%     res(j+1) = H(j+1,j)*abs(y(j))/norm(b);
    
%     j = j+1;
end


% y = H(1:m,1:m)\(beta*e1);
x = x0 + V(:,1:j)*y;


return