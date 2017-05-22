function [ FS,CLASS ] = featureExtractionSSVEP( X, range, filtRange, cont, start, Fs )
%featureExtractionSSVEP
    % Input Vars:
%  X            = Signal Vector;
%  range        = size of windows to test (e.g. 250 ? 2500 in 60 point increments)
%  filtRange    = Bandpass range for butterworth filter
%  cont         = Plot Graph if Empty, else run to completion
%  start        = Data point to start analyzing dataset from
%  Fs           = Sampling Rate (250 Hz by default)

    % Default Values:
if nargin < 6
    Fs = 250;
end
if nargin < 5
    start = 1;
end
if nargin < 4
    cont = 0;
end
if nargin < 3
    filtRange = [8,20];
end
if nargin < 2
    range = 250:60:2500;
end

    %Preallocate:
% coder.varsize('F');
% F = [];
NUMBERFEATURES = 40;
F  = zeros(size(range,2),NUMBERFEATURES);
FS = zeros(1,NUMBERFEATURES*size(range,2));
P = zeros(size(range,2),56);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fprintf('Current index = [%d to %d]\r\n',start, fin);
    fprintf('length = %d\r\n',range(i));
    fch = customFilt(X(start:fin),Fs,filtRange,3);
        %%%Feature Extraction: (per channel)
    [F(i,:),P(i,:)] = fESSVEP(fch,Fs,isempty(cont));
%     P(i,:) = fESSVEP2(fch,Fs,isempty(cont));
    if isempty(cont)
        commandwindow;
        cont = input('Approve/continue?\n');
        clf(12);
    end
end
idx = 1:4;
M = zeros(1,4);
L = M;

if ~isempty(cont)
    f = P(1,1:28);
    figure(13);hold on;xlim([8 20]);
    for i = 1:4
        [M(i),L(i)] = max(max(P(:,29+((i-1)*7):35+((i-1)*7))));
        if(cont~=0)
            plot(f((i-1)*7+L(i)),M(i),'or');
        end
    end
    [Peak,ClusterLoc] = max(M);
    idx2 = idx(ClusterLoc~=idx);
    b = zeros(1,length(idx2));
    for i=1:length(idx2)
        b(i) = (M(idx2(i)) > Peak/3);
    end
    h = refline([0,Peak/3]); h.Color = 'r';
    for i = 1:size(P,1)
        plot(P(i,1:28),P(i,29:end))
    end
    commandwindow;
    if(cont~=0)
        CLASS = input('Approve/continue?\n');
    else
        CLASS = [];
    end
    if isempty(CLASS)
        if sum(b)==0
            CLASS = ClusterLoc
        else
            CLASS = 0
        end
    end
    clf(13);
else
    CLASS = 0;
end

for i = 1:size(F,1)
    FS(1+size(F,2)*(i-1):size(F,2)*(i)) = F(i,:);
end

% FS=F(:);
end

