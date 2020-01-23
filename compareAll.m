function allresults = compareAll()

%Creating matrix A from the dataset : 'Abalone RBF kernel with sigma = .15'
% A is a PSD matrix
load abalone_dataset

input.sigma = .15;% this is the sparsity parameter, small sigam means we want sprse matrix, sigma=1 give a dense matrix.
input.datasetbasename = 'Abalone';
input.k = 20; %k, the desired rank of the approximation
input.chunk = 10; % % how often to reorthogonalize in Spectral and  Power-method

Ab = CreateDistMatrix(abaloneInputs');
input.A = zeros(size(Ab));
for row=1:size(Ab,1)
    input.A(row, :) = exp(-Ab(row,:)/input.sigma^2);
end
clear Ab;
%sampNumbers is a vector of the numbers of column samples to use
% We are going to run for 20 differnt number of sampling:
sampNumbers =[20    27    34    41    48    55    62    68    75    82    89    96   103   109   116   123   130   137   144   150];%
input.sampleNum = sampNumbers; 
allresults.input = input;

%% Calling diffrent implementation of Nystrom method here. First we comput the exact leverage scores and optimal rank-k approximation error.
tic
[U, Sigma] = sparseSorteig(input.A, input.k);%get the first k+1 sorted eigenvalues.
U1t = U(:, 1:input.k)';
allresults.levscores = sum(U1t.*U1t);
input.exactlevscoretiming = toc;
input.levscoreprobs = allresults.levscores/input.k;
topspectrum = diag(Sigma(1:input.k,1:input.k));
%The optimal error (by Frobenius norm):
allresults.err = sqrt(norm(input.A, 'fro')^2 - sum(topspectrum.^2));
%Calling all of the Nystrom methods for each one of the different number of
%smaplings
for sn = 1:length(input.sampleNum)
    input.l = input.sampleNum(sn);    
    unif_Data(sn) = NystromUniform(input);
    exac_Data(sn) = NystromExactLev(input);
    spec_Data(sn) = NystromSpectralLev(input);
    powe_Data(sn) = NystromPowerLev(input);
    frob_Data(sn) = NystromFrobLev(input);   
end
allresults.unif_Data = unif_Data;
allresults.exac_Data = exac_Data;
allresults.spec_Data = spec_Data;
allresults.powe_Data = powe_Data;
allresults.frob_Data = frob_Data;

save('allresults');
clear in;


