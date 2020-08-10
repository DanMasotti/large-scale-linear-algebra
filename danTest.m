clear all; close all; clc;

disp('--- Vanilla GMRES -----');
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');

guess = zeros(rows,1);
m = 80;

tic
[x,y]= myGMRES(A,b,guess,80,b);
toc

figure
fig = semilogy(y/norm(b,2));
ylim([10^-15 1])
xlabel('$$Iterations$$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$$Residual$$', 'interpreter', 'latex', 'fontsize', 16);
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 16);
grid on;
title('GMRES on Sherman Matrix', 'interpreter', 'latex', 'fontsize', 20);
grid on;

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
[x,y]= myGMRES(A_,b_,guess,m,b);
toc

figure
fig = semilogy(y);
xlabel('$$Iterations$$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$$Residual$$', 'interpreter', 'latex', 'fontsize', 16);
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 16);
grid on;
title('GMRES on 10-diagonal Sherman Matrix', 'interpreter', 'latex', 'fontsize', 20);
err_my_GMRES = norm(b_-A_*x,2)

 

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
[x,y]= myGMRES(A_,b_,guess,m,b);
toc

figure
fig = semilogy(y);
xlabel('$$Iterations$$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$$Residual$$', 'interpreter', 'latex', 'fontsize', 16);
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 16);
grid on;
title('GMRES on 4-diagonal Sherman Matrix', 'interpreter', 'latex', 'fontsize', 20);
ylim([10^-14 10^2]);
err_my_GMRES = norm(b_-A_*x,2)

