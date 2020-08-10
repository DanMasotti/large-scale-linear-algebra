clear all; close all; clc;

%WEEK OF 3/4 - Lanczos, CG Method

% n = 1000;
% %sparse symmetric matrix
% A = zeros(n,n);
% A = A + 4*eye(n);
% offset = zeros(n-1,n-1) + eye(n-1,n-1);
% A(2:n,1:n-1) = A(2:n,1:n-1) + offset;
% A(1:n-1,2:n) = A(1:n-1,2:n) + offset;
% guess = zeros(n,1);
% b = A*ones(n,1);
% 
% 
% disp('----- Lanczos Iteration -----')
% m = 10;
% v1 = ones(n,1);
% [T, V] = my_lanczos(A,v1,m);
% %TODO: find a way to measure error
% error = norm(T-(V'*A*V))
% 
% disp('---- Lanczos for Linear Systems ----')
% m = 10;
% [y,x]= lanczos_solver(A,guess,b,m);
% figure
% semilogy(y)
% title('Lanczos: Residuals VS. Iterations');
% err = norm(b-A*x,2)
%  
% 
% 
% %FOR WEEK 3/10 -- FINISH CG, Ch 10
% 
disp('---- Conjugate Gradient ----')
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');
m = 10;
guess = zeros(rows,1);
[y,x] = my_cg(A,guess,b,m);
figure
semilogy(y)
title('CG: Residuals VS. Iterations');
err = norm(b-A*x,2)

%spy(A)**

%FOR WEEK 3/17 -- Implement SSOR, GS, ILU, preconditioner (Ch 10) and try on
%Sherman (GMRES), update proposal to update Background (intuition about
%inverse matrix being sum)
disp('---- Comparison of GMRES on Sherman with Preconditioning ----');
A = mmread('sherman2.mtx');
n = size(A,1);
b = mmread('sherman2_rhs1.mtx');
% disp('-- none --');
% tic;
% x = gmres(A,b);
% toc
% err = norm(x-A\b,2)
% % [x , A\b]
% 
% % disp('-- ILU --');
% % [L,U,R] = my_ilu(A);
% % tic;
% % x = gmres(L*U,b);
% % toc
% % disp('me')
% % err = norm(x-A\b,2)
% % % [x , A\b]
% % 
% % disp('matlabs ILU')
% % [L, U] = ilu(A);
% % x = gmres(L*U,b);
% % disp('them')
% % err = norm(x-A\b,2)
% % % [x A\b]
% 
% disp('-- GS --');
% % d = diag(A);
% D = diag(diag(A));
% E = -tril(A,-1);
% F = -triu(A,1);
% M = (D-E)*(eye(n)\D)*(D-F);
% A_ = (eye(n)\M)*A;
% b_ = (eye(n)\M)*b;
% tic
% x = gmres(A_,b_,n,1e-12,n);
% toc
% err = norm(x-A\b,2)
% % [x, A\b]

% disp('-- SSOR --');
% D = diag(diag(A));
% E = -tril(A,-1);
% F = -triu(A,1);
% w = 1/1024;
% M = (D-w*E)*(inv(D))*(D-w*F)';
% A_ = (eye(n)\M)*A;
% b_ = (eye(n)\M)*b;
% tic
% x = gmres(A,b, [],1e-15,500,M);
% toc
% err = norm(x-A\b,2)
% [x , A\b]

%FOR WEEK 4/4, fix the code, use the GRMES on the linear solve
%(using matlab's iterative methods)

%FOR WEEK 4/11, time the time stepper of back-slash vs gmres.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,rows,cols] = mmread('sherman2.mtx');
b = mmread('sherman2_rhs1.mtx');
% view structure of A
figure
spy(A);

% let's construct a preconditioner from the first 10 diagonals of A
d = -5:1:5;
B = spdiags(A,d); 
%grabs the 10 diagonals of A 
M = spdiags(B,d,rows,cols); 
%construct a new matrix of the 10 diagonals of A
Af = full(A);
e = eig(Af); % compute eigenvalues
figure
fig = plot(real(e),imag(e),'r*');
xlabel('$$Re(\lambda)$$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$$Im(\lambda)$$', 'interpreter', 'latex', 'fontsize', 16);
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 16);
grid on;
title('Spectrum of $$A$$ before Preconditioning', 'interpreter', 'latex', 'fontsize', 20);
saveas(fig,strcat('before','.png'))


A = M\A;
b = M\b;

% convert A to full matrix format
Af = full(A);
e = eig(Af); % compute eigenvalues
figure
fig = plot(real(e),imag(e),'r*');
xlabel('$$Re(\lambda)$$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$$Im(\lambda)$$', 'interpreter', 'latex', 'fontsize', 16);
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 16);
grid on;
title('Spectrum of $$A$$ with Preconditioner', 'interpreter', 'latex', 'fontsize', 20);
saveas(fig,strcat('after','.png'))
