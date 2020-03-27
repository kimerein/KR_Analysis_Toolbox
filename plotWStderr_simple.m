function [data1,data2,isBigEnough]=plotWStderr_simple(data1,data2,trialDuration,c1,c2,times,allLines)

norm=1;
% normWindow=[4 4.7];
% normWindow=[3 6];
normWindow=[1 3];
baseSubtract=1;
baseWindow=[0 0.65];
% baseWindow=[0 2.3];
% baseWindow=[0.5 1.2];
% baseWindow=[0 3];
% baseWindow=[0.3 1.25];
% baseWindow=[0.5 0.78];
% baseWindow=[1.5 3];
% baseWindow=[3.7 4]; % kim's expts
% baseWindow=[2.385 2.685]; % arbora's expts
% baseWindow=[2.7 2.9]; % arbora's expts
% baseWindow=[2 2.5]; % arbora's expts
% baseWindow=[1 2]; % arbora's expts
ds=1;
doFill=1;
isBigEnough=[];

t=times;

if baseSubtract==1
    for i=1:size(data1,1)
        base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        data1(i,:)=data1(i,:)-base1;
        data2(i,:)=data2(i,:)-base2;
    end
end

% thresh=2*10^4; 
% thresh=0.6*10^4; 
% thresh=0.25*10^4;
thresh=0.05*10^4; 
% thresh=3.5*10^4;
% thresh=0.44*10^4; % for 4 Hz
% thresh=0.10*10^4; % for 14 Hz
% thresh=0.6*10^4; % for 20 Hz
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh) & nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>0 & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>0;
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% isBigEnough=nanmean(data1(:,t>=5 & t<=6),2)>thresh;
isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh & nanmean(data2(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
% isBigEnough=(nanmean(data1(:,t>=normWindow(1) & t<=normWindow(2)),2)>thresh);
data1=data1(isBigEnough,:);
data2=data2(isBigEnough,:);
disp(sum(isBigEnough==1));
 
baseWindow=[0.51 1]; % for LOW FREQ
baseWindow=[0.9 1]; % for HIGH FREQ

if norm==1
    for i=1:size(data1,1)
        scale=nanmean(data1(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data1(i,:)=data1(i,:)./scale;
%         scale=nanmean(data2(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data2(i,:)=data2(i,:)./scale;
    end
end

if ~isempty(ds)
    data1=downSampMatrix(data1,ds);
    data2=downSampMatrix(data2,ds);
    if length(trialDuration)>1
        trialDuration=downSampAv(trialDuration,ds);
    else
        trialDuration=downSampAv(t,ds);
    end
end

if length(trialDuration)>1 & allLines~=1
    dropthebase=[0 1];
%     dropthebase=[3 3.7];
%     dropthebase=[0 4];
%     dropthebase=[0 3];
%     average_baseline_data1=nanmean(nanmean(data1(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
%     average_baseline_data2=nanmean(nanmean(data2(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
    average_baseline_data1=0;
    average_baseline_data2=0;
    figure(); 
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    hold on;
    if doFill==1
        fill([trialDuration fliplr(trialDuration)],[nanmean(data1,1)-average_baseline_data1+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmean(data1,1)-average_baseline_data1-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
    end
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2,'Color',c2);
    if doFill==1
        fill([trialDuration fliplr(trialDuration)],[nanmean(data2,1)-average_baseline_data2+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmean(data2,1)-average_baseline_data2-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
    end
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2,'Color',c2);
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    plot(trialDuration,nanmean(data2,1)-average_baseline_data2-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    
    plot(trialDuration,nanmean(data1,1)-average_baseline_data1,'Color',c1);
    
elseif length(trialDuration)>1 & allLines==1
    
    figure(); 
    plot(trialDuration,data1','Color',c1);
    hold on;
    plot(trialDuration,data2','Color',c2);
    
else
    figure();
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    hold on;
    if doFill==1
        fill([linspace(0,trialDuration,size(data1,2)) fliplr(linspace(0,trialDuration,size(data1,2)))],[nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
    end
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
    
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
    if doFill==1
        fill([linspace(0,trialDuration,size(data2,2)) fliplr(linspace(0,trialDuration,size(data2,2)))],[nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
    end
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1),'Color',c2);
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    plot(linspace(0,trialDuration,size(data2,2)),nanmean(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
    
    plot(linspace(0,trialDuration,size(data1,2)),nanmean(data1,1),'Color',c1);
    
end
end