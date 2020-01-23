function result = NystromUniform(input)
%This method implements Nystrom with uniform sampling.
%  -input.A,  PSD matrix
%  -input.k, the desired rank of the approximation
%  -input.l , the number of columns to sample and
n = size(input.A,1);
tic
perm = randperm(n);
perm = perm(1:input.l);
C = input.A(:, perm);
W = full(C(perm,:));
[V,D] = sorteig(W);%eigen decomposition
Wkinv = V(:, 1:input.k)*pinv(D(1:input.k, 1:input.k))*V(:,1:input.k)';
result.timings= toc;
result.err = norm(input.A - C*Wkinv*C', 'fro');
end