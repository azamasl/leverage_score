function result = NystromExactLev(input)
% This method implements Nystrom with exact sketch leverage scores
%  -input.A,  PSD matrix
%  -input.k, the desired rank of the approximation
tic
n = size(input.A,1);
Id = speye(n);
colindices = ones(1,input.l);
%fprintf('This is levscoreprobs from exact');
%mydataE.levscoreprobs =input.levscoreprobs;
for i=1:input.l
    colindices(i) = find(cumsum(input.levscoreprobs) >= rand(),1);%picks elements of levscoreprob randomly
end

%fprintf('This is colindices ');
%mydataE.colindices = colindices;
scaling = input.levscoreprobs(colindices).^(1/2);%computing the scaling factors
%mydataE.scaling=scaling;
S = Id(:,colindices)*diag(1./scaling);
%fprintf('This is S ');
%mydataE.S = S;

C = input.A*S;

W = full(S'*C);
%mydata.W = W;
%save('mydataE');
[V,D] = sorteig(W);
Wkinv = V(:, 1:input.k)*pinv(D(1:input.k, 1:input.k))*V(:,1:input.k)';
result.timings = input.exactlevscoretiming + toc;
result.err = norm(input.A - C*Wkinv*C', 'fro');
end
