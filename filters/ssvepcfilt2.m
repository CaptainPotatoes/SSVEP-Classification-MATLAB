function [ Y ] = ssvepcfilt2( X )
% [5 40] bandpass butterworth N=3
X=X(:);
b = [0.0418768282347742,0,-0.125630484704323,0,0.125630484704323,0,-0.0418768282347742];
a = [1,-3.99412602172993,6.79713743558926,-6.44840721730666,3.65712515526032,-1.17053739881085,0.159769122451512];
Y = filtfilt(b,a,X);
end
