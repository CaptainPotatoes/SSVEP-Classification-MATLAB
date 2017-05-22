function [ FS,CLASS ] = classifySSVEP( X, start, Fs )
%CLASSIFYSSVEP - FINAL VERSION FOR MATLAB CODER
    % INPUT VARS:
    % X - input array (any size)
    % start - where to start from in 'X'
    % Fs - signal sampling frequency
    
% range - range of window sizes to view
range = 250:60:1000; % 1-4 s at 60pt intervals
NUMFTS = 40;
NUMP = 56;
F  = zeros(size(range,2),NUMFTS);
FS = zeros(1,NUMFTS*size(range,2));
P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
%     fprintf('Current index = [%d to %d]\r\n',start, fin);
%     fprintf('length = %d\r\n',range(i));
%     fch = customFilt(X(start:fin),Fs,filtRange,3);
    fch = ssvepcfilt(X(start:fin));
        %%%Feature Extraction: (per channel)
    [F(i,:),P(i,:)] = fESSVEP(fch,Fs,false);
end
idx = 1:4;
M = zeros(1,4);
L = M;
% f = P(1,1:28);
% figure(13);hold on;xlim([8 20]);
for i = 1:4
    [M(i),L(i)] = max(max(P(:,29+((i-1)*7):35+((i-1)*7))));
%     plot(f((i-1)*7+L(i)),M(i),'or');
end
[Peak,ClusterLoc] = max(M);
idx2 = idx(ClusterLoc~=idx);
b = zeros(1,length(idx2));
for i=1:length(idx2)
    b(i) = (M(idx2(i)) > Peak/3);
end

for i = 1:size(P,1)
%     plot(P(i,1:28),P(i,29:end))
end

if sum(b)==0
    CLASS = ClusterLoc;
else
    CLASS = 0;
end

for i = 1:size(F,1)
    FS(1+size(F,2)*(i-1):size(F,2)*(i)) = F(i,:);
end

end

