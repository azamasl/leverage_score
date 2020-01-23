function [V,D] = sorteig(M, varargin)
% Returns sorted eigenvalues
[V,D] = eig(M, varargin{:});
[~,idx] = sort(diag(D),1,'descend');
V = V(:, idx);
D = diag(D);
D = diag(D(idx));

end