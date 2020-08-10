function [T, V] = my_lanczos(A,v_1, m)
threshold = 10^-6;
%preallocation
n = size(A,1);
alphas = zeros(m,1);
betas = zeros(m,1);
V = zeros(n,m);
T = eye(m);
V(:,1) = v_1/norm(v_1);

%step j = 1
w = A*V(:,1);
alphas(1) = w'*V(:,1);
w = w - alphas(1)*V(:,1);
betas(1)=0;
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
end
V = V(:,1:m);
T = T(1:m,1:m);
end