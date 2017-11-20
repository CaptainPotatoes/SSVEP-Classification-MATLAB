function yfit = knn(tsX, tX, tY, Knn)
%function yfit = knnclassification(testsamplesX,samplesX, samplesY, Knn, type)
% Classify using the Nearest neighbor algorithm
% Inputs:
% 	tX	   - Train samples
%	tY	   - Train labels
%   tsX (testsamplesX) - Test  samples to classify
%	Knn		       - Number of nearest neighbors 
%
% Outputs
%	result	- Predicted targets
%if nargin < 5
%    type = '2norm';
%end

L			= length(tY);
Uc          = unique(tY);

if (L < Knn)
   error('You specified more neighbors than there are points.')
end

N                 = size(tsX, 1);
yfit              = zeros(N,1);

for i = 1:N
    dist            = sum((tX - ones(L,1)*tsX(i,:)).^2,2);
    [~, indices]    = sort(dist);  
    n               = hist(tY(indices(1:Knn)), Uc);
    [~, best]       = max(n);
    yfit(i)         = Uc(best);
end

end
