function [levscores, it] = power_method_approx_levscores(input)
%This methods approximates the lev scorese with power method.
%  -A,  PSD matrix
%  -k, the desired rank of the approximation
n = size(input.A,1);
l = input.k;
S = randn(n, l);% S is a matrix of nxl of random numbers
C = input.A*S;
it = 0;
tol = 0.01;
levscores = zeros(1,n);
while it < input.maxiters
    it = it + 1;
    if (rem(it,input.chunk)==0)%if we need to reorthogonalize
        [Q,~] = qr(C, 0);
        C = Q;
        prev_levscores = levscores;
        levscores = sum(Q.^2,2)';% sum of the squares for each column of Q
        if ( norm(levscores - prev_levscores, Inf) < tol)
            return;
        end
    end
    C = input.A*(input.A*C); % since A is PSD, C = A*A^T*C
end
