function [ sigs ] = generateTestSignal( f_new, len )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
sigs = zeros(length(f_new),len);
for i = 1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),len);
end
sigs = sigs(:)';

end

