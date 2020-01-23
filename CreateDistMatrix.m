function distance_matrix = CreateDistMatrix(datapoints)
%Each xi represents a datapoint and correseponds to a row in 'datapoints'
% This method returns the distance matrix, that is: ||xi-xj||_2^2 as the (ij)th entry
% First: center the data and makes std dev of points 1. This is called
% whiten the noise.
numpoints = size(datapoints,1);
m = mean(datapoints);
stdv = std(datapoints);
stdv(stdv < 10^(-4)) = 1; % avoid division by 0
pts = (datapoints - repmat(m, numpoints, 1))./repmat(stdv, numpoints, 1);
distance_matrix = zeros(numpoints, numpoints);
parfor row = 1:numpoints % do for in parallel
    displacements = (pts - repmat(pts(row, :),numpoints,1))';
    distance_matrix(row, :) = sum(displacements.^2);
end

