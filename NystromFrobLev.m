function result = NystromFrobLev(input)
% This method implements Nystrom with approximating frobenius sketch leverage scores
%  -input.A,  PSD matrix
%  -input.k, the desired rank of the approximation
n = size(input.A,1);
Id = speye(n);
tic
r = 2*input.k;
S = randn(n, r); %Gaussian sampling
[Q,~] = qr(input.A*S, 0); % produces unitary matrix Q and UT matrix R so that A*S = Q*R.
mydataF.Q = Q;
[U,~,~] = svds(Q'*input.A, input.k);
levscores = sum((Q*U).^2, 2)';%sums along the columns so levscores is a vector.
mydataF.levscores= levscores;
levscoreprob = levscores/input.k;% a vector
% sample according to the computed leverage scores above:
indices = ones(1,input.l);% Column indices.
for i=1:input.l
    indices(i) = find(cumsum(levscoreprob) >= rand(),1);% returns the index of the first levscoreprob it find which is larger than the rand()
end
mydataF.indices = indices;
scaling = levscoreprob(indices).^(1/2);%computing the scaling factors
mydataF.scaling = scaling;
S = Id(:,indices)*diag(1./scaling);% From the sparse identity matrix of size nxn pics the columns with indices and scale them. So S is nxl.
mydataF.S = S;
C = input.A*S;
W = full(S'*C); %W is lxl .
mydataF.W = W;
save('mydataF');
[V,D] = sorteig(W);%eigen decomposition
Wkinv = V(:, 1:input.k)*pinv(D(1:input.k, 1:input.k))*V(:,1:input.k)';
result.timings = toc;
result.err = norm(input.A - C*Wkinv*C', 'fro');
end

