clear all; close all; clc;

% disp('-------SPARSE TEST----------');
% A = zeros(100,100);
% A = A + 4*eye(100);
% offset = zeros(99,99) - eye(99,99);
% A(2:100,1:99) = A(2:100,1:99) + offset;
% A(1:99,2:100) = A(1:99,2:100) + offset;
% guess = zeros(100,1);
% b = A*ones(100,1);
% disp('--- Vinilla GMRES -----');
% [y,x]= my_GMRES(A,b,guess,20);
% figure
% semilogy(y)
% title('Sparse Test with GMRES')
% err_my_GMRES = norm(b-A*x,2)
% disp('---- Restarted my_GMRES -----');
% [y,x]= restarted_my_GMRES(A,b,guess,20);
% figure
% semilogy(y)
% title('Sparse Test with GMRES Restarted')
% err_restart = norm(b-A*x,2)
% 
% disp('------RANDOM TEST----------');
% A = rand(100,100);
% b = rand(100,1);
% guess = rand(100,1);
% m = 100;
% disp('--- Vinilla GMRES -----');
% [y,x] = my_GMRES(A,b,guess,m);
% figure
% semilogy(y)
% title('Random Test with GMRES')
% title('GMRES residuals')
% err_my_GMRES = norm(b-A*x,2)
% disp('---- Restarted my_GMRES -----');
% [y,x]= restarted_my_GMRES(A,b,guess,20);
% figure
% semilogy(y)
% title('Random Test with GMRES Restarted')
% err_restart = norm(b-A*x,2)

disp('-------SHERMAN TEST---------');
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');

% view structure of A
% figure
% spy(A);

% let's construct a preconditioner from the first 10 diagonals of A
d = -5:1:5;
B = spdiags(A,d); %grabs the 10 diagonals of A 
M = spdiags(B,d,rows,cols); %construct a new matrix of the 10 diagonals of A


A_ = M\A;
b_ = M\b;

% convert A to full matrix format
% Af = full(A);
% e = eig(Af); % compute eigenvalues
% figure
% plot(real(e),imag(e),'r*');

guess = zeros(rows,1);

m = 80;
% m_res = 7;
disp('--- Vanilla GMRES -----');
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');
tic
[y,x]= my_GMRES(A,b,guess,m);
toc
figure
fig = semilogy(y);
xlabel('Iterations');
ylabel('Residual');
grid on;
title('GMRES on Sherman Matrix')

err_my_GMRES = norm(b-A*x,2)

disp('--- 10-diagGMRES -----');
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');
% let's construct a preconditioner from the first 10 diagonals of A
d = -5:1:5;
B = spdiags(A,d); %grabs the 10 diagonals of A 
M = spdiags(B,d,rows,cols); %construct a new matrix of the 10 diagonals of A
A_ = M\A;
b_ = M\b;
tic
[y,x]= my_GMRES(A_,b_,guess,m);
toc
figure
fig = semilogy(y);
xlabel('Iterations');
ylabel('Residual');
grid on;
title('GMRES on 10-diagonal Sherman Matrix')
err_my_GMRES = norm(b-A*x,2)

disp('---------4-diagGMRES-----------')
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');
d = -2:1:2;
n = size(A,1);
B = spdiags(A,d);
M = spdiags(B,d,n,n);
A_ = M\A;
b_ = M\b;

tic
[y,x]= my_GMRES(A_,b_,guess,m);
toc
figure
fig = semilogy(y);
xlabel('Iterations');
ylabel('Residual');
grid on;
title('GMRES on 4-diagonal Sherman Matrix')
err_my_GMRES = norm(b-A*x,2)

disp('---- Restarted my_GMRES -----');
%preconditioning
% d= -6:1:6;
% n = size(A,1);
% B = spdiags(A,d);
% M = spdiags(B,d,n,n);
% A = M\A;
% b = M\b;
% tic
% [y,x]= restarted_my_GMRES(A_,b_,guess,m_res);
% toc
% figure
% semilogy(y)
% title('Restarted GMRES on Preconditioned Sherman')
% err_restart = norm(b-A*x,2)

function [Y,x] = my_GMRES(A,b,guess,m)

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
    I = eye(j+1);
    H_tilda = H(1:j+1,1:j);
    %y = ls_helper(H_tilda,beta*I(:,1));
    y = H_tilda\(beta*I(:,1));
    r_j = (norm(beta*I(:,1)-H_tilda*y))/norm(b);
    Y(j) = r_j;
    if r_j <= 10^-16
        m = j;
        break;
    end
    if (abs(H(j+1,j)) < 10^-6)
        disp('breaks down'); 
        x = x_0 + V(:,1:j)*y;
        return
    end
    V(:,j+1)=w_j/H(j+1,j);
end
x = x_0 + V(:,1:m)*y;
end

function [Y,x] = restarted_my_GMRES(A,b,guess,m)
% m = int32(m/2);

Y = [];

xi = guess;
for batch = 1:10
    [Yi,xi] = my_GMRES(A,b,xi,m);
    
    Y = [Y; Yi];
end

x = xi;
end
