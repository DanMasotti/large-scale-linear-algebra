function [N,x] = my_cg(A,x_0,b,m)
n = size(A,1);
R = zeros(n,m);
x = x_0;
R(:,1) = b - A*x;
p = R(:,1);
for j = 1:m
    a = (R(:,j)'*R(:,j))/((A*p)'*p);
    x = x + a*p;
    R(:,j+1) = R(:,j)-a*A*p;
    beta = (R(:,j+1)'*R(:,j+1))/(R(:,j)'*R(:,j));
    p = R(:,j+1) + beta*p;
end
N = zeros(1,m);
for i = 1:m
    N(1,i) = norm(R(:,i),2);
end
end 