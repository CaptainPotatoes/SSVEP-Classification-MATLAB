function yfit = TestMasterNew(thetanull)
%
% Run input data computation ----------------------------------------------

%arduin = serial('COM7','Baudrate',38400);
%fopen(arduin);

%Variables
s = 10000; %Initial matrix size (total data size)
%WindowMin = 9900; %XLim min
%WindowMax = 10000; %XLim max
%thetanull=zeros(s,6); %Don't uncomment this
%thetanull(1:end-1)=thetanull(2:end);
theta2null=zeros(s,6); %Don't uncomment this
persistent humanActivityDataTestMinusActivity
persistent humanActivityDataTestMinusActivity2
%persistent Activity
persistent ActivityModified
XAccSample = zeros(500,20);
YAccSample = zeros(500,20);
ZAccSample = zeros(500,20);
AccResSample = zeros(500,20);
XGyroSample = zeros(500,20);
YGyroSample = zeros(500,20);
ZGyroSample = zeros(500,20);
GyroResSample = zeros(500,20);

%
if isempty(humanActivityDataTestMinusActivity2)


    theta2Test=thetaTest(:,1:3).^2;
    theta3Test=theta2Test(:,1)+theta2Test(:,2)+theta2Test(:,3);
    theta4Test=theta3Test.^(1/2);
    theta5Test=thetaTest(:,4:6).^2;
    theta6Test=theta5Test(:,1)+theta5Test(:,2)+theta5Test(:,3);
    theta7Test=theta6Test.^(1/2);

    XAccTest = zeros(500,20);
    YAccTest = zeros(500,20);
    ZAccTest = zeros(500,20);
    AccResTest = zeros(500,20);
    XGyroTest = zeros(500,20);
    YGyroTest = zeros(500,20);
    ZGyroTest = zeros(500,20);
    GyroResTest = zeros(500,20);

     for k = 1:500
        for p = 1 + 20*(k-1)
            for j = 20*k
                XAccTest(k,:) = thetaTest(p:j,1);
                YAccTest(k,:) = thetaTest(p:j,2);
                ZAccTest(k,:) = thetaTest(p:j,3);
                AccResTest(k,:) = theta4Test(p:j);
                XGyroTest(k,:) = thetaTest(p:j,4);
                YGyroTest(k,:) = thetaTest(p:j,5);
                ZGyroTest(k,:) = thetaTest(p:j,6);
                GyroResTest(k,:) = theta7Test(p:j);
            end
        end
     end
%{
    Activity = zeros(500,1);
    Activity(1:6,1) = 1;
    Activity(7,1) = 3;
    Activity(8,1) = 1 ;
    Activity(9,1) = 5;
    Activity(10:25,1) = 1 ;
    Activity(26:27,1) = 3 ;
    Activity(28:43,1) = 1;
    Activity(44,1) = 3 ;
    Activity(45,1) = 5 ;
    Activity(46:60,1) = 1; 
    Activity(61,1) = 5 ;
    Activity(62:104,1) = 1; 
    Activity(105,1) = 3; 
    Activity(106:126,1) = 1; 
    Activity(127,1) = 5 ;
    Activity(128:141,1) = 1 ;
    Activity(142,1) = 5 ;
    Activity(143:157,1) = 1 ;
    Activity(158,1) = 5 ;
    Activity(159:211,1) = 1 ;
    Activity(212,1) = 5 ;
    Activity(213:222,1) = 1 ;
    Activity(223,1) = 5;
    Activity(224:234,1) = 1 ;
    Activity(235,1) = 5 ;
    Activity(236:243,1) = 1 ;
    Activity(244,1) = 5;
    Activity(245:255,1) = 1 ;
    Activity(256,1) = 5 ;
    Activity(257:285,1) = 1 ;
    Activity(286,1) = 5 ;
    Activity(287:291,1) = 1 ;
    Activity(292,1) = 5 ;
    Activity(293:298,1) = 1 ;
    Activity(299,1) = 5  ;
    Activity(300:306,1) = 1 ;
    Activity(307,1) = 5 ;
    Activity(308:312,1) = 1 ;
    Activity(313,1) = 5 ;
    Activity(314:321,1) = 1 ;
    Activity(322:393,1) = 3 ;
    Activity(394:436,1) = 4;
    Activity(437:467,1) = 2;
    Activity(468:479,1) = 6;
    Activity(480:485,1) = 1;
    Activity(486:500,1) = 2;
%}
    ActivityModified=zeros(32,1);
    ActivityModified(1,1)=5;
    ActivityModified(2,1)= 5 ;
    ActivityModified(3,1)= 5 ;
    ActivityModified(4,1) = 5 ;
    ActivityModified(5,1)= 5 ;
    ActivityModified(6,1) = 5 ;
    ActivityModified(7,1) = 5 ;
    ActivityModified(8,1) = 5;
    ActivityModified(9,1) = 5 ;
    ActivityModified(10,1) = 5;
    ActivityModified(11,1) = 5 ;
    ActivityModified(12,1) = 5 ;
    ActivityModified(13,1) = 5 ;
    ActivityModified(14,1) = 5  ;
    ActivityModified(15,1) = 5 ;
    ActivityModified(16,1) = 5 ;
    ActivityModified(17:24,1) = 1 ;
    ActivityModified(25:32,1)= 1 ;
    
    T_mean = horzcat(Wmean(XAccTest),Wmean(YAccTest), Wmean(ZAccTest), Wmean(AccResTest), Wmean(XGyroTest), Wmean(YGyroTest), Wmean(ZGyroTest), Wmean(GyroResTest));
    T_stdv = horzcat(Wstd(XAccTest),Wstd(YAccTest), Wstd(ZAccTest), Wstd(AccResTest), Wstd(XGyroTest), Wstd(YGyroTest), Wstd(ZGyroTest), Wstd(GyroResTest));
    T_pca  = horzcat(Wpca1(XAccTest),Wpca1(YAccTest), Wpca1(ZAccTest), Wpca1(AccResTest), Wpca1(XGyroTest), Wpca1(YGyroTest), Wpca1(ZGyroTest), Wpca1(GyroResTest));
    T_max = horzcat(Wmax(XAccTest),Wmax(YAccTest), Wmax(ZAccTest), Wmax(AccResTest), Wmax(XGyroTest), Wmax(YGyroTest), Wmax(ZGyroTest), Wmax(GyroResTest));
    T_min = horzcat(Wmin(XAccTest),Wmin(YAccTest), Wmin(ZAccTest), Wmin(AccResTest), Wmin(XGyroTest), Wmin(YGyroTest), Wmin(ZGyroTest), Wmin(GyroResTest));
    T_countmin = horzcat(WCountMin(XAccTest),WCountMin(YAccTest), WCountMin(ZAccTest), WCountMin(AccResTest), WCountMin(XGyroTest), WCountMin(YGyroTest), WCountMin(ZGyroTest), WCountMin(GyroResTest));
    T_countmax = horzcat(WCountMax(XAccTest),WCountMax(YAccTest), WCountMax(ZAccTest), WCountMax(AccResTest), WCountMax(XGyroTest), WCountMax(YGyroTest), WCountMax(ZGyroTest), WCountMax(GyroResTest));
    T_Integrate = horzcat(WIntegrate(XAccTest),WIntegrate(YAccTest), WIntegrate(ZAccTest), WIntegrate(AccResTest), WIntegrate(XGyroTest), WIntegrate(YGyroTest), WIntegrate(ZGyroTest), WIntegrate(GyroResTest));

    %humanActivityDataTest = horzcat(T_mean, T_stdv, T_pca, T_max, T_min, T_countmin, T_countmax, T_Integrate, Activity);
    humanActivityDataTestMinusActivity = horzcat(T_mean, T_stdv, T_pca, T_max, T_min, T_countmin, T_countmax, T_Integrate);
    humanActivityDataTestMinusActivity2 = [humanActivityDataTestMinusActivity(9,:);humanActivityDataTestMinusActivity(45,:);humanActivityDataTestMinusActivity(61,:);humanActivityDataTestMinusActivity(127,:);humanActivityDataTestMinusActivity(142,:);humanActivityDataTestMinusActivity(158,:);humanActivityDataTestMinusActivity(212,:);humanActivityDataTestMinusActivity(223,:);humanActivityDataTestMinusActivity(235,:);humanActivityDataTestMinusActivity(244,:);humanActivityDataTestMinusActivity(256,:);humanActivityDataTestMinusActivity(286,:);humanActivityDataTestMinusActivity(292,:);humanActivityDataTestMinusActivity(299,:);humanActivityDataTestMinusActivity(307,:);humanActivityDataTestMinusActivity(313,:);humanActivityDataTestMinusActivity(314:321,:);humanActivityDataTestMinusActivity(322:329,:)];
else
end
%end
%
%uicontrol('string','stop','callback','jin=0;');
delay=0.01;
%jin=1;

%while(jin)
    %thetanull(end,:)=fscanf(arduin,'%d %d %d %d %d %d');
    
    theta2null(:,1:3) = thetanull(:,1:3)./16384;
    theta2null(:,4:6) = thetanull(:,4:6)./131;
    theta = theta2null;
     
    theta2=theta(:,1:3).^2;
    theta3=theta2(:,1)+theta2(:,2)+theta2(:,3);
    theta4=theta3.^(1/2);
    
    theta5=theta(:,4:6).^2;
    theta6=theta5(:,1)+theta5(:,2)+theta5(:,3);
    theta7=theta6.^(1/2);
    
    thetanull(1:end-1)=thetanull(2:end);
    theta2null(1:end-1)=theta2null(2:end);
    theta(1:end-1)=theta(2:end);
    
    for k = 1:500
        for p = 1 + 20*(k-1)
            for j = 20*k
                XAccSample(k,:) = theta(p:j,1);
                YAccSample(k,:) = theta(p:j,2);
                ZAccSample(k,:) = theta(p:j,3);
                AccResSample(k,:) = theta4(p:j);
                XGyroSample(k,:) = theta(p:j,4);
                YGyroSample(k,:) = theta(p:j,5);
                ZGyroSample(k,:) = theta(p:j,6);
                GyroResSample(k,:) = theta7(p:j);
            end
        end
    end

    %plotGraph=plot(theta4);
    %set(plotGraph,'YData',theta4)
    %ylim([-5 5]);
    %xlim([WindowMin WindowMax]);
    %time=linspace(0,2,100);
    %xlim([time(1) time(end)]);

if max(max(theta4(9900:end)))>2
    %Run humanactivitylearnertest-----------------------------------------

    %1 = Idle (Sitting or Standing)
    %2 = Transitioning
    %3 = Walk
    %4 = Run
    %5 = Fall
    %6 = Jumping
    %{
    N=500;
    for z = 1:N
        for q = 1+20*(N-1)
            for r = 20*N
                if max(max(theta4(q:r)))>2
    %}
    T_meanSample = horzcat(Wmean(XAccSample),Wmean(YAccSample), Wmean(ZAccSample), Wmean(AccResSample), Wmean(XGyroSample), Wmean(YGyroSample), Wmean(ZGyroSample), Wmean(GyroResSample));
    T_stdvSample = horzcat(Wstd(XAccSample),Wstd(YAccSample), Wstd(ZAccSample), Wstd(AccResSample), Wstd(XGyroSample), Wstd(YGyroSample), Wstd(ZGyroSample), Wstd(GyroResSample));
    T_pcaSample  = horzcat(Wpca1(XAccSample),Wpca1(YAccSample), Wpca1(ZAccSample), Wpca1(AccResSample), Wpca1(XGyroSample), Wpca1(YGyroSample), Wpca1(ZGyroSample), Wpca1(GyroResSample));
    T_maxSample = horzcat(Wmax(XAccSample),Wmax(YAccSample), Wmax(ZAccSample), Wmax(AccResSample), Wmax(XGyroSample), Wmax(YGyroSample), Wmax(ZGyroSample), Wmax(GyroResSample));
    T_minSample = horzcat(Wmin(XAccSample),Wmin(YAccSample), Wmin(ZAccSample), Wmin(AccResSample), Wmin(XGyroSample), Wmin(YGyroSample), Wmin(ZGyroSample), Wmin(GyroResSample));
    T_countminSample = horzcat(WCountMin(XAccSample),WCountMin(YAccSample), WCountMin(ZAccSample), WCountMin(AccResSample), WCountMin(XGyroSample), WCountMin(YGyroSample), WCountMin(ZGyroSample), WCountMin(GyroResSample));
    T_countmaxSample = horzcat(WCountMax(XAccSample),WCountMax(YAccSample), WCountMax(ZAccSample), WCountMax(AccResSample), WCountMax(XGyroSample), WCountMax(YGyroSample), WCountMax(ZGyroSample), WCountMax(GyroResSample));
    T_IntegrateSample = horzcat(WIntegrate(XAccSample),WIntegrate(YAccSample), WIntegrate(ZAccSample), WIntegrate(AccResSample), WIntegrate(XGyroSample), WIntegrate(YGyroSample), WIntegrate(ZGyroSample), WIntegrate(GyroResSample));

    humanActivityDataSample = horzcat(T_meanSample, T_stdvSample, T_pcaSample, T_maxSample, T_minSample, T_countminSample, T_countmaxSample, T_IntegrateSample);

%Run TestMaster--------------------------------------------------------
%yfit = trainedClassifierNew3.predictFcn(humanActivityDataSample);
%yfit(1:end-1)=yfit(2:end);
%yfit=zeros(500,1);
%

    yfit = knnclassification(humanActivityDataSample(end-9:end,:),humanActivityDataTestMinusActivity2,ActivityModified,3);
    if max(max(yfit)) == 5 && sum(yfit(end-4:end,1))==5
        alpha=1
    end
else
    yfit = zeros(10,1);
    %----------------------------------------------------------------------
    %{
    if max(max(yfit)) == 5  && sum(yfit(496:end)) == 5
           Status = 'Fall';
           txt = uicontrol('string',Status,...
           'Position', [750 50 600 150],'FontWeight','Bold','BackgroundColor','red','FontSize',20);
    end
    %}
%else
end
pause(delay);
end
%end

%fclose(arduin);