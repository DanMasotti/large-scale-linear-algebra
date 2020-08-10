A = rand(100,100)*100+100;
b = rand(100,1)*100 + 100;
iters = 100;
N = 100;
D = [zeros(1,N/2),ones(1,N/2)];
D = D(randperm(2*N/2))';
A = randn(N,N);
P = A * diag(D)*inv(A);
[hess,kSpace]  = arnoldi(P,b,iters);

err = norm(kSpace(:,1:end-1)'*P*kSpace(:,1:end-1)-hess(1:end-1,:), 2)


function [H,V] = arnoldi(A,b,m)
n = size(A,1);
v_i = b/norm(b);
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = v_i;
for j = 1:m
    for i = 1:j
       H(i,j) = (A*V(:,j))'*V(:,i);
    end
    
    w_j = A*V(:,j);
    for i=1:j
        w_j = w_j - H(i,j)*V(:,i);
    end
    H(j+1,j)=norm(w_j,2);
   if (abs(H(j+1,j)) < 10^-6)
      disp('breaks down');
      return
   end
   V(:,j+1)=w_j/H(j+1,j);
   
end

end