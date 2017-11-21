function [ Accuracy ] = runFileClassification( flnm, figs,range,START )
[DATA, ~] = csvread(flnm);
if(isempty(figs))
    figs = false;
end
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
range_cca = 250;
maxlen = max(range);
winHop = 60;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==10
THRESHOLD_FRACTION = 2 ; i = 1;
% Test Signals:
sti_f = [ 9.5, 15.15, 16.67, 18.52, 20.0 ];
t_length=4;                              % data length (4 s)
TW=1:1:t_length;
TW_p=round(TW*Fs);
N=2;    % number of harmonics
ref0=refsig(sti_f(1),Fs,t_length*Fs,N);
ref1=refsig(sti_f(2),Fs,t_length*Fs,N);
ref2=refsig(sti_f(3),Fs,t_length*Fs,N);
ref3=refsig(sti_f(4),Fs,t_length*Fs,N);
ref4=refsig(sti_f(5),Fs,t_length*Fs,N);
PredictedClass = zeros(length(range),length(wStart));
PC2 = PredictedClass; PC3 = PC2; ActualClass = PC3; 
% %{ 
% frequencies = (0:(998/2-1))*250/998;
for i = 1:length(wStart)
    i
    for j = 1:length(range)
        fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(range(j)-1));
        [Y(i,:), Ppsd(i,:,j), featurePoint(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)), PLOTDATA, THRESHOLD_FRACTION, START);
%         [Y(i,:), PSDfull(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)), THRESHOLD_FRACTION, PLOTDATA);
%         [P, F] = extractPowerSpectrum2ch(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)));
%         if(PLOTDATA)
%             figure(1); 
%             plot(f,P(i,:));
%         end
        if(size(DATA,2) == 3)
            ActualClass(j,i) = mode(DATA(wStart(i):wStart(i)+(range(j)-1),3));
        end
        X12 = [X1(wStart(i):wStart(i)+(range_cca(j)-1)), X2(wStart(i):wStart(i)+(range_cca(j)-1))]';
        tw_length = range_cca(j)/Fs;
        [wx0,wy0,r0]=cca(X12,ref0(:,1:TW_p(tw_length)));
        [wx1,wy1,r1]=cca(X12,ref1(:,1:TW_p(tw_length)));
        [wx2,wy2,r2]=cca(X12,ref2(:,1:TW_p(tw_length)));
        [wx3,wy3,r3]=cca(X12,ref3(:,1:TW_p(tw_length)));
        [wx4,wy4,r4]=cca(X12,ref4(:,1:TW_p(tw_length)));
        [v,idx]=max([max(r0),max(r1),max(r2),max(r3),max(r4)]);
        PC3(j,i) = idx-1;
    end
end
Accuracy = zeros(1,length(range)); Accuracy2 = Accuracy; Accuracy3 = Accuracy;

if(size(DATA,2) == 3)
    for j = 1:length(range)
        if(figs)
            figure(j); hold on;
            plot(ActualClass(j,:));
            plot(Y(:,1),'.');
            plot(Y(:,2),'*');
            plot(PC3(j,:),'.');
            legend('Actual Class','Y1','Y2','CCA PC3'); xlabel('Classifier Calls (1/s)');
            ylabel('Class Label');
%         title(['SSVEP Classifier Output (', num2str(range(j)) ,' samples) ; Accuracy: ' num2str(100*max(Accuracy(:,j))) '%']);
            title(['SSVEP Classifier Output (', num2str(range(j)) ,' samples) ; Accuracy: ' num2str(100*Accuracy(1)) '%']);
        end
        Compare1 = ActualClass(j,:) == Y(:,1)';% | ActualClass == PredictedClass2;
        Compare2 = ActualClass(j,:) == Y(:,2)';
        Compare3 = ActualClass(j,:) == PC3(j,:);
        Accuracy(1,j) = sum(Compare1)/length(ActualClass(j,:));
        Accuracy(2,j) = sum(Compare2)/length(ActualClass(j,:));
        Accuracy(3,j) = sum(Compare3)/length(ActualClass(j,:));
        
    end
else
    if(figs)
        plot(PredictedClass); plot(PC2);
        legend('Predicted Class', 'PC2'); xlabel('Classifier Calls (1/s)');
        ylabel('Class Label');
        title(['SSVEP Classifier Output; Accuracy: ??%']);
    end
end
Accuracy



end

