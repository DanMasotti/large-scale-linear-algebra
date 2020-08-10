function [R, x] = lanczos_solver(A,x_0,b,m)
threshold = 10^-6;
r = b - A*x_0;
beta = norm(r,2);

%preallocation
I = eye(m);
n = size(A,1);
alphas = zeros(m,1);
betas = zeros(m,1);
V = zeros(n,m);
T = zeros(m);
V(:,1) = r/norm(r);
R = zeros(m);

%step j = 1
w = A*V(:,1);
alphas(1) = w'*V(:,1);
w = w - alphas(1)*V(:,1);
% betas(1)=r;
T(1,1) = alphas(1);

%cycle
for j = 1:m
   if j == 1
       w = A*V(:,j);
   else 
      w = A*V(:,j)-betas(j)*V(:,j-1);
   end  
   alphas(j) = w'*V(:,j);
   w = w-alphas(j)*V(:,j);
   betas(j+1) = norm(w,2);
   if betas(j+1) <= threshold
       disp('breaks down');
       m = j;
       break;
   end
   V(:,j+1) = w/betas(j+1);
   T(j,j) = alphas(j);
   T(j+1,j) = betas(j+1);
   T(j,j+1) = betas(j+1);
   y = T(1:j,1:j)\(beta*I(1:j,1));
   R(j) = betas(j+1)*norm(I(1:j,j)'*y);
end
V = V(:,1:m);
T = T(1:m,1:m);
x = x_0 + V(:,1:m)*y;
end