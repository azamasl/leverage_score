function result = NystromPowerLev(input)
tic
n = size(input.A,1);
result.num1s = 0;
input.maxiters = ceil(log(1/3)/log(.9)); % Assuming the spectral gap ratio is smaller than .9
Id = speye(n);
% compute the approximate leverage scores with power method
[result.approxlevscores, result.num1s] = power_method_approx_levscores(input);
levscoreprobs = result.approxlevscores/input.k;
% sample according to those leverage scores
colindices = ones(1,input.l);
for i=1:input.l
    colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
end
scaling = levscoreprobs(colindices).^(1/2);%computing the scaling factors
S = Id(:,colindices)*diag(1./scaling);
C = input.A*S;
W = full(S'*C);
[V,D] = sorteig(W);
Wkinv = V(:, 1:input.k)*pinv(D(1:input.k, 1:input.k))*V(:,1:input.k)';
result.timings= toc;
result.err = norm(input.A - C*Wkinv*C', 'fro');
end
