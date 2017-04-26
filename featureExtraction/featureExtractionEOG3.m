function [ F ] = featureExtractionEOG3( X, LTH1, LTH2, UTH1, UTH2, plotData )
%% FeatureExtractionEOG3 (from diff signal)
% Accepts 2 low threshold and 2 upper thresholds.
% Upper Thresholds UTH1 UTH2
% Lower Threshold. LTH1 LTH2
%%-TODO: FINALIZE THRESHOLDS:
% UTH1 = 0.4E-4;
% UTH2 = 2.75E-4;
% LTH1 = -0.5E-4;
% LTH2 = -2.75E-4;
[Fmax, Imax] = max(X);
[Fmin, Imin] = min(X);
Famplitude = Fmax-Fmin;
Fstd = std(X);
FInt1 = trapz(X(X>UTH1 & X<UTH2));
FInt2 = trapz(-X(X>LTH2 & X<LTH1));
Fvelocity = Famplitude/(Imax - Imin);
FcountMin = sum(X>LTH2 & X<LTH1); %Count between bottom two lines.
FcountMax = sum(X>UTH1 & X<UTH2); 
FcountMaxHigh = sum(X>UTH2);
FcountMinLow = sum(X<LTH2);
%%% PLOT FEATURES %%%
if plotData
    figure(4);hold on;
    IDX = 1:250;
    plot(IDX(X>LTH2 & X<LTH1),X(X>LTH2 & X<LTH1),'k.');
    plot(IDX(X>UTH1 & X<UTH2),X(X>UTH1 & X<UTH2),'k^');
    plot(Imax,Fmax,'r*');
    plot(Imin,Fmin,'r*');
end

F = [Fmax,Fmin,Famplitude,Fstd,FInt1,FInt2,Fvelocity,FcountMin,FcountMax,FcountMaxHigh,FcountMinLow];

end

%{ 
EOF
%}
