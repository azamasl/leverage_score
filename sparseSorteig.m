function [V,D] = sparseSorteig(M, varargin)

[V,D] = eigs(M, varargin{:});
[~,idx] = sort(diag(D),1,'descend');
V = V(:, idx);
D = diag(D);
D = diag(D(idx));

end