clear all;clc;
arduin = serial('COM13','Baudrate',115200);
fopen(arduin);
%s=10000;
thetanull=zeros(100,6);
uicontrol('string','stop','callback','jin=0;');
delay=0.01;
jin=1;
yfit = zeros(5,1);
while(jin)
%     thetanull(end,:)=fscanf(arduin,'%d %d %d %d %d %d');
    theta = fscanf(arduin, '%d %d')
    %{
    thetanull(1:end-1)=thetanull(2:end);
    [yfit] = TestMasterNew(thetanull);   
    theta2null(:,1:3) = thetanull(:,1:3)./16384;
    %theta2null(:,4:6) = thetanull(:,4:6)./131;
    %theta = theta2null;
     
    %theta2=theta(:,1:3).^2;
    %theta3=theta2(:,1)+theta2(:,2)+theta2(:,3);
    %theta4=theta3.^(1/2);
    
    %theta5=theta(:,4:6).^2;
    %theta6=theta5(:,1)+theta5(:,2)+theta5(:,3);
    %theta7=theta6.^(1/2);
     for k=1:100
        theta4(k,1)=norm(theta2null(k,1:3));
        %theta7(k,1)=norm(theta2null(k,4:6));
     end
    
    plot(theta4);
    ylim([-0 5]);
    %xlim([9900 10000]);
    yfit
    %alpha;
    %}
end
fclose(arduin);