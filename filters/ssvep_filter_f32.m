% Optimized SSVEP C Filter:

function [ Y ] = ssvep_filter_f32( X, Fs )
% [5 40] bandpass butterworth N=3; 250Hz
% X = X(:);
Y = single(X)
if (Fs == 250)
b = [0.0418768282347742,0,-0.125630484704323,0,0.125630484704323,0,-0.0418768282347742];
a = [1,-3.99412602172993,6.79713743558926,-6.44840721730666,3.65712515526032,-1.17053739881085,0.159769122451512];
Y = filtfilt(b,a,X);
Y = single(Y);
elseif Fs == 500

elseif Fs == 1000
end


end

