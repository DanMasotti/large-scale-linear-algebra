clear all; close all; clc;

disp('-------SPARSE TEST----------');
A = zeros(100,100);
A = A + 4*eye(100);
offset = zeros(99,99) - eye(99,99);
A(2:100,1:99) = A(2:100,1:99) + offset;
A(1:99,2:100) = A(1:99,2:100) + offset;
guess = zeros(100,1);
b = A*ones(100,1);
disp('--- Vinilla FOM -----');
[y,x]= fom(A,b,guess,20);
figure
semilogy(y)
title('Sparse Test with FOM');
err_fom = norm(b-A*x,2)
disp('---- Restarted FOM -----');
[y,x]= restarted_fom(A,b,guess,20);
figure
semilogy(y)
title('Sparse Test with Restarted FOM');
err_restart = norm(b-A*x,2)

disp('------RANDOM TEST----------');
A = rand(100,100);
b = rand(100,1);
guess = rand(100,1);
m = 100;
disp('--- Vinilla FOM -----');
[y,x] = fom(A,b,guess,m);
figure
semilogy(y)
title('Random Test with FOM');
err_fom = norm(b-A*x,2)
disp('---- Restarted FOM -----');
[y,x]= restarted_fom(A,b,guess,20);
figure
semilogy(y)
title('Random Test with Restarted FOM');
err_restart = norm(b-A*x,2)

disp('-------SHERMAN TEST---------');
[A,rows,cols] = mmread('sherman2.mtx');
guess = zeros(rows,1);
b = mmread('sherman2_rhs1.mtx');
m = cols;
disp('--- Vinilla FOM -----');
[y,x]= fom(A,b,guess,800);
figure
semilogy(y)
title('Sherman Test with FOM');
err_fom = norm(b-A*x,2)
disp('---- Restarted FOM -----');
[y,x]= restarted_fom(A,b,guess,20);
figure
semilogy(y)
title('Sherman Test with Restarted FOM');
err_restart = norm(b-A*x,2)

function [Y,x] = fom(A,b,guess,m)
n = size(A,1);
x_0 = guess;
r_0 = b - A*x_0;
beta = norm(r_0,2);
v_i = r_0/beta;
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = v_i;
Y = zeros(m);
for j = 1:m
    w_j = A*V(:,j);
    for i = 1:j
        H(i,j) = w_j'*V(:,i);
        w_j = w_j - H(i,j)*V(:,i);
    end  
    H(j+1,j) = norm(w_j,2);
    I = eye(j);
    y = H(1:j,1:j)\(beta*I(:,1));
    r_j = (H(j+1,j)*abs(y(j)))/norm(b);
    Y(j) = r_j; 
    if (abs(H(j+1,j)) < 10^-6)
        disp('breaks down'); 
        x = x_0 + V(:,1:j)*y;
        return
    end
    V(:,j+1)=w_j/H(j+1,j);
end
x = x_0 + V(:,1:m)*y;
end

function [Y, x] = restarted_fom(A,b,guess,m)
m = int32(m/4);
for batch = 1:100
    [Y,x] = fom(A,b,guess,m);
end
end
