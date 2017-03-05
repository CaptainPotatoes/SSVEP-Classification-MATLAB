function yfit = knn(testsamplesX,samplesX, samplesY, Knn)
%function yfit = knnclassification(testsamplesX,samplesX, samplesY, Knn, type)
% Classify using the Nearest neighbor algorithm
% Inputs:
% 	samplesX	   - Train samples
%	samplesY	   - Train labels
%   testsamplesX   - Test  samples
%	Knn		       - Number of nearest neighbors 
%
% Outputs
%	result	- Predicted targets
%if nargin < 5
%    type = '2norm';
%end

L			= length(samplesY);
Uc          = unique(samplesY);

if (L < Knn),
   error('You specified more neighbors than there are points.')
end

N                 = size(testsamplesX, 1);
yfit              = zeros(N,1);

for i = 1:N,
    dist            = sum((samplesX - ones(L,1)*testsamplesX(i,:)).^2,2);
    [m, indices]    = sort(dist);  
    n               = hist(samplesY(indices(1:Knn)), Uc);
    [m, best]       = max(n);
    yfit(i)         = Uc(best);
end

end
