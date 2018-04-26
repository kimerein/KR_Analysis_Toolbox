function output=plotF1analysis(datadir,trialDuration)

if iscell(datadir)
    all_noTheta_trialAv_noLED=[];
    all_noTheta_trialAv_LED=[];
    all_theta_trialAv_noLED=[];
    all_theta_trialAv_LED=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'noTheta_trialAv_noLED']);
        noTheta_trialAv_noLED=a.noTheta_trialAv;
        all_noTheta_trialAv_noLED=concatStructs(all_noTheta_trialAv_noLED,noTheta_trialAv_noLED);
        
        a=load([d '\' 'noTheta_trialAv_LED']);
        noTheta_trialAv_LED=a.noTheta_trialAv;
        all_noTheta_trialAv_LED=concatStructs(all_noTheta_trialAv_LED,noTheta_trialAv_LED);
        
        a=load([d '\' 'theta_trialAv_noLED']);
        theta_trialAv_noLED=a.theta_trialAv;
        all_theta_trialAv_noLED=concatStructs(all_theta_trialAv_noLED,theta_trialAv_noLED);

        a=load([d '\' 'theta_trialAv_LED']);
        theta_trialAv_LED=a.theta_trialAv;
        all_theta_trialAv_LED=concatStructs(all_theta_trialAv_LED,theta_trialAv_LED);
    end
    noTheta_trialAv_noLED=all_noTheta_trialAv_noLED;
    noTheta_trialAv_LED=all_noTheta_trialAv_LED;
    theta_trialAv_noLED=all_theta_trialAv_noLED;
    theta_trialAv_LED=all_theta_trialAv_LED;
else 
    a=load([datadir '\' 'noTheta_trialAv_noLED']);
    noTheta_trialAv_noLED=a.noTheta_trialAv;
    
    a=load([datadir '\' 'noTheta_trialAv_LED']);
    noTheta_trialAv_LED=a.noTheta_trialAv;
    
    a=load([datadir '\' 'theta_trialAv_noLED']);
    theta_trialAv_noLED=a.theta_trialAv;
    
    a=load([datadir '\' 'theta_trialAv_LED']);
    theta_trialAv_LED=a.theta_trialAv;
end

% noTheta_trialAv_noLED.F1amp=noTheta_trialAv_noLED.F1amp(1:92,:)+noTheta_trialAv_noLED.F1amp(93:184,:)+noTheta_trialAv_noLED.F1amp(185:276,:)+noTheta_trialAv_noLED.F1amp(277:368,:);
% noTheta_trialAv_noLED.F1amp=noTheta_trialAv_noLED.F1amp/4;
% theta_trialAv_noLED.F1amp=theta_trialAv_noLED.F1amp(1:92,:)+theta_trialAv_noLED.F1amp(93:184,:)+theta_trialAv_noLED.F1amp(185:276,:)+theta_trialAv_noLED.F1amp(277:368,:);
% theta_trialAv_noLED.F1amp=theta_trialAv_noLED.F1amp/4;
% noTheta_trialAv_LED.F1amp=noTheta_trialAv_LED.F1amp(1:92,:)+noTheta_trialAv_LED.F1amp(93:184,:)+noTheta_trialAv_LED.F1amp(185:276,:)+noTheta_trialAv_LED.F1amp(277:368,:);
% noTheta_trialAv_LED.F1amp=noTheta_trialAv_LED.F1amp/4;
% theta_trialAv_LED.F1amp=theta_trialAv_LED.F1amp(1:92,:)+theta_trialAv_LED.F1amp(93:184,:)+theta_trialAv_LED.F1amp(185:276,:)+theta_trialAv_LED.F1amp(277:368,:);
% theta_trialAv_LED.F1amp=theta_trialAv_LED.F1amp/4;

% noTheta_trialAv_noLED.F1amp=noTheta_trialAv_noLED.F1amp(1:92,:)+noTheta_trialAv_noLED.F1amp(93:184,:)+noTheta_trialAv_noLED.F1amp(185:276,:);
% noTheta_trialAv_noLED.F1amp=noTheta_trialAv_noLED.F1amp/3;
% theta_trialAv_noLED.F1amp=theta_trialAv_noLED.F1amp(1:92,:)+theta_trialAv_noLED.F1amp(93:184,:)+theta_trialAv_noLED.F1amp(185:276,:);
% theta_trialAv_noLED.F1amp=theta_trialAv_noLED.F1amp/3;
% noTheta_trialAv_LED.F1amp=noTheta_trialAv_LED.F1amp(1:92,:)+noTheta_trialAv_LED.F1amp(93:184,:)+noTheta_trialAv_LED.F1amp(185:276,:);
% noTheta_trialAv_LED.F1amp=noTheta_trialAv_LED.F1amp/3;
% theta_trialAv_LED.F1amp=theta_trialAv_LED.F1amp(1:92,:)+theta_trialAv_LED.F1amp(93:184,:)+theta_trialAv_LED.F1amp(185:276,:);
% theta_trialAv_LED.F1amp=theta_trialAv_LED.F1amp/3;

[data1,data2]=plotWStderr(noTheta_trialAv_noLED.F1amp(1:end,:),theta_trialAv_noLED.F1amp(1:end,:),trialDuration,'b','r');
t=linspace(0.5,trialDuration-0.5,size(data1,2));
[data1,data2]=plotWStderr(noTheta_trialAv_noLED.F1amp(1:end,:),noTheta_trialAv_LED.F1amp(1:end,:),trialDuration,'b','c');
% t=linspace(0.165,trialDuration-0.165,size(data1,2));
% t=linspace(0.5,trialDuration-0.5,size(data1,2));
temp1=nanmean(data1(:,t>1.25 & t<=1.35),2);
temp2=nanmean(data2(:,t>1.25 & t<=1.35),2);
output.noTheta_early=(nanmean(temp2)-nanmean(temp1))./nanmean(temp1);
% isBigEnough=temp1>nanmedian(temp1) | temp2>nanmedian(temp2);
% temp1=temp1(isBigEnough==1);
% temp2=temp2(isBigEnough==1);
% disp(nanmean(temp1))
% output.noTheta_con=nanmean(temp1);
% disp(nanstd(temp1,[],1)./sqrt(length(temp1)))
% output.noTheta_con_se=nanstd(temp1,[],1)./sqrt(length(temp1));
% disp(nanmean(temp2))
% output.noTheta_led=nanmean(temp2);
% disp(nanstd(temp2,[],1)./sqrt(length(temp2)))
% output.noTheta_led_se=nanstd(temp2,[],1)./sqrt(length(temp2));
% output.noTheta_p=signrank(temp1,temp2);
% disp(output.noTheta_p)
% output.noTheta_LEDminusCon=temp2-temp1;
% output.noTheta_LEDminusCon_divideByCon=(temp2-temp1)./temp1;
temp1=nanmean(data2(:,t>2.4 & t<=2.5),2);
temp2=nanmean(data1(:,t>2.4 & t<=2.5),2);
output.noTheta_late=(nanmean(temp2)-nanmean(temp1))./nanmean(temp1);
% % output.noTheta_early=(temp1-temp2)./(temp1+temp2);
% output.noTheta_early=(nanmean(temp1)-nanmean(temp2))/nanmean(temp2);
% % temp1(temp1<0)=0;
% % temp2(temp2<0)=0;
% % output.noTheta_early(isinf(output.noTheta_early))=nan;

% temp1=nanmean(data2(:,t>2.3 & t<=2.8),2);
% temp2=nanmean(data1(:,t>2.3 & t<=2.8),2);
% output.noTheta_late=(nanmean(temp1)-nanmean(temp2))/nanmean(temp2);
% % temp1(temp1<0)=0;
% % temp2(temp2<0)=0;
% output.noTheta_late(isinf(output.noTheta_late))=nan;

[data1,data2]=plotWStderr(theta_trialAv_noLED.F1amp(1:end,:),theta_trialAv_LED.F1amp(1:end,:),trialDuration,'r','c');
temp1=nanmean(data1(:,t>1.1 & t<=1.2),2);
temp2=nanmean(data2(:,t>1.1 & t<=1.2),2);
output.theta_early=(nanmean(temp1)-nanmean(temp2))./nanmean(temp1);
% isBigEnough=temp1>nanmedian(temp1) | temp2>nanmedian(temp2);
% temp1=temp1(isBigEnough==1);
% temp2=temp2(isBigEnough==1);
% disp(nanmean(temp1))
% output.theta_con=nanmean(temp1);
% disp(nanstd(temp1,[],1)./sqrt(length(temp1)))
% output.theta_con_se=nanstd(temp1,[],1)./sqrt(length(temp1));
% disp(nanmean(temp2))
% output.theta_led=nanmean(temp2);
% disp(nanstd(temp2,[],1)./sqrt(length(temp2)))
% output.theta_led_se=nanstd(temp2,[],1)./sqrt(length(temp2));
% output.theta_p=signrank(temp1,temp2);
% disp(output.theta_p)
% output.theta_LEDminusCon=temp2-temp1;
% output.theta_LEDminusCon_divideByCon=(temp2-temp1)./temp1;
% temp1=nanmean(data2(:,t>1.2 & t<=1.7),2);
% temp2=nanmean(data1(:,t>1.2 & t<=1.7),2);
% output.theta_early=(nanmean(temp2)-nanmean(temp1))/nanmean(temp2);
% % temp1(temp1<0)=0;
% % temp2(temp2<0)=0;
% % output.theta_early(isinf(output.theta_early))=nan;
% 
temp1=nanmean(data2(:,t>2.4 & t<=2.5),2);
temp2=nanmean(data1(:,t>2.4 & t<=2.5),2);
output.theta_late=(nanmean(temp1)-nanmean(temp2))./nanmean(temp1);
% output.theta_late=(nanmean(temp2)-nanmean(temp1))/nanmean(temp2);
% % temp1(temp1<0)=0;
% % temp2(temp2<0)=0;
% % output.theta_late(isinf(output.theta_late))=nan;

% [data1,data2]=plotWStderr(noTheta_trialAv_LED.F1amp(1:end,:),theta_trialAv_LED.F1amp(1:end,:),trialDuration,'r','c'); 
% [data1,data2]=plotWStderr(noTheta_trialAv_noLED.F1amp,theta_trialAv_LED.F1amp,trialDuration,'r','c');
 
output.noTheta_trialAv_noLED.F1amp=noTheta_trialAv_noLED.F1amp; 
output.theta_trialAv_noLED.F1amp=theta_trialAv_noLED.F1amp;
output.noTheta_trialAv_LED.F1amp=noTheta_trialAv_LED.F1amp;
output.theta_trialAv_LED.F1amp=theta_trialAv_LED.F1amp;

disp('snoop dogg');
return

if iscell(datadir)
    all_noTheta_noLED=[];
    all_noTheta_LED=[];
    all_theta_noLED=[];
    all_theta_LED=[];
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'noTheta_noLED']);
        noTheta_noLED=a.noTheta;
        all_noTheta_noLED=concatStructs(all_noTheta_noLED,noTheta_noLED);
        
        a=load([d '\' 'noTheta_LED']);
        noTheta_LED=a.noTheta;
        all_noTheta_LED=concatStructs(all_noTheta_LED,noTheta_LED);
        
        a=load([d '\' 'theta_noLED']);
        theta_noLED=a.theta;
        all_theta_noLED=concatStructs(all_theta_noLED,theta_noLED);

        a=load([d '\' 'theta_LED']);
        theta_LED=a.theta;
        all_theta_LED=concatStructs(all_theta_LED,theta_LED);
    end
    noTheta_noLED=all_noTheta_noLED;
    noTheta_LED=all_noTheta_LED;
    theta_noLED=all_theta_noLED;
    theta_LED=all_theta_LED;
else
    a=load([datadir '\' 'noTheta_noLED']);
    noTheta_noLED=a.noTheta;
    
    a=load([datadir '\' 'noTheta_LED']);
    noTheta_LED=a.noTheta;
    
    a=load([datadir '\' 'theta_noLED']);
    theta_noLED=a.theta;
    
    a=load([datadir '\' 'theta_LED']);
    theta_LED=a.theta;
end

[data1,data2]=plotWStderr(noTheta_noLED.F1amp,theta_noLED.F1amp,trialDuration,'b','r');
[noTheta_noLED.F1amp,noTheta_LED.F1amp]=plotWStderr(noTheta_noLED.F1amp,noTheta_LED.F1amp,trialDuration,'b','c');
[data3,data4]=plotWStderr(theta_noLED.F1amp,theta_LED.F1amp,trialDuration,'r','c');
% plotWStderr([noTheta_noLED.F1amp; theta_noLED.F1amp],[noTheta_noLED.F1amp; theta_noLED.F1amp],trialDuration,'k','k');

output.noTheta_noLED.F1amp=noTheta_noLED.F1amp;
output.theta_noLED.F1amp=theta_noLED.F1amp;
output.noTheta_LED.F1amp=noTheta_LED.F1amp;
output.theta_LED.F1amp=theta_LED.F1amp;

disp('sir mixalot');

end

function newStruct=concatStructs(struct1,struct2)

if isempty(struct1)
    newStruct=struct2;
    return
elseif isempty(struct2)
    newStruct=struct1;
    return
end
fie=fieldnames(struct1);
for i=1:length(fie)
    f=fie{i};
    if isfield(struct2,f)
        newStruct.(f)=[struct1.(f); struct2.(f)];
    end
end

end

% function plotWStderr(data1,data2,trialDuration,c1,c2)
% 
% figure(); 
% plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
% hold on;
% % fill([linspace(0,trialDuration,size(data1,2)) fliplr(linspace(0,trialDuration,size(data1,2)))],[nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
% plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
% plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
% plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
% 
% plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
% % fill([linspace(0,trialDuration,size(data2,2)) fliplr(linspace(0,trialDuration,size(data2,2)))],[nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
% plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
% plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
% plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
% end