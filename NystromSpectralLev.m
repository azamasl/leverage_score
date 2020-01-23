function result = NystromSpectralLev(input)
tic
n = size(input.A,1);
Id = speye(n);
S = randn(n, 2*input.k);
q = 4;
Y = input.A*S;
iter = 0; % keep track of how many iterations of the power method were used
while iter < q
    iter = iter + 1;    
    if (rem(iter,input.chunk)==0)% reorthogonalize every chunk steps
        [Q,~] = qr(Y, 0);
        Y = Q;
    end
    Y = input.A*(input.A*Y); % since A is PSD, Y = A*A^T*Y
end

[U,~,~] = svds(Y, input.k);
result.approxlevscores = sum(U.^2, 2)';
levscoreprobs = result.approxlevscores/input.k;
% sample according to those leverage scores
colindices = ones(1,input.l);
for i=1:input.l
    colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
end
scalingfactors = levscoreprobs(colindices).^(1/2);
S = Id(:,colindices)*diag(1./scalingfactors);
C = input.A*S;
W = full(S'*C);
[V,D] = sorteig(W);
Wkinv = V(:, 1:input.k)*pinv(D(1:input.k, 1:input.k))*V(:,1:input.k)';
result.timings = toc;
result.err = norm(input.A - C*Wkinv*C', 'fro');
end
