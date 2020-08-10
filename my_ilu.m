function [L,U,R] = my_ilu(A)
threshold = 10^-6;
n = size(A,1);
original = A;

for i = 1:n
    for k = 1:i-1
        if A(i,k) > threshold
            A(i,k) = A(i,k)/A(k,k);
            for j = k+1:n
                if A(i,j) > threshold
                    A(i,j) = A(i,j) - A(i,k)*A(k,j);
                end
            end
        end
    end
end
L = tril(A);
U = triu(A);
R = L*U-original;
end