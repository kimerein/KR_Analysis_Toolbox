function compareBestAndWorstResponses(psth,downSamp,responseFromPSTH,ampOrder_max,ampOrder_min)

% assumes all units from same experiment, recorded simultaneously

t=psth.t;
stimcond=psth.unitStimcond{1};
ledcond=psth.unitLED{1};

u=unique(stimcond);

for i=1:length(u)
    currstimcond=u(i);
    [black,red]=averageUnitPSTHs(responseFromPSTH,u(i));
    psthByStimcond.(['black' num2str(i)])=black;
    psthByStimcond.(['red' num2str(i)])=red;
end

for i=1:length(u)
    temp=psthByStimcond.(['black' num2str(i)]);
    avResponseByStimcond.(['black' num2str(i)])=nanmean(temp(:,t>=4 & t<=6.5),2);
    temp=psthByStimcond.(['red' num2str(i)]);
    avResponseByStimcond.(['red' num2str(i)])=nanmean(temp(:,t>=4 & t<=6.5),2);
end

maxInd=nan(1,size(psthByStimcond.black1,1));
minInd=nan(1,size(psthByStimcond.black1,1));
for i=1:size(psthByStimcond.black1,1)
    % by units
    responses=[];
    for j=1:length(u)
        % by stimcond
        temp=avResponseByStimcond.(['black' num2str(j)]);
        responses=[responses temp(i)];
    end
    % find max response
    [~,maxInd(i)]=nanmax(responses);
    % find min response
    [~,minInd(i)]=nanmin(responses);
end

if ~isempty(ampOrder_max)
    maxInd=ampOrder_max;
    minInd=ampOrder_min;    
end

for i=1:length(u)
    currstimcond=u(i);
    [black,red]=averageUnitPSTHs(psth,u(i));
    psthByStimcond_forFigure.(['black' num2str(i)])=black;
    psthByStimcond_forFigure.(['red' num2str(i)])=red;
end

% make best and worst response psths
for i=1:size(psthByStimcond_forFigure.black1,1)
    % by units
    temp=psthByStimcond_forFigure.(['black' num2str(maxInd(i))]);
    bestPSTH_black(i,:)=temp(i,:);
    temp=psthByStimcond_forFigure.(['black' num2str(minInd(i))]);
    worstPSTH_black(i,:)=temp(i,:);
    
    temp=psthByStimcond_forFigure.(['red' num2str(maxInd(i))]);
    bestPSTH_red(i,:)=temp(i,:);
    temp=psthByStimcond_forFigure.(['red' num2str(minInd(i))]);
    worstPSTH_red(i,:)=temp(i,:);
end

downt=downSampAv(t,downSamp);

sub_plotWStderr(downSampMatrix(bestPSTH_black,downSamp),downSampMatrix(worstPSTH_black,downSamp),downt,'c','k');
title('No LED best vs worst');

sub_plotWStderr(downSampMatrix(bestPSTH_red,downSamp),downSampMatrix(worstPSTH_red,downSamp),downt,'r','m');
title('LED best vs worst');


end


function [data1,data2]=sub_plotWStderr(data1,data2,trialDuration,c1,c2)

norm=0;
% normWindow=[4.5 6];
normWindow=[2.8 3];
normWindow=[0 3];
baseSubtract=0;
% baseWindow=[0 2.5];
baseWindow=[0 4];
% baseWindow=[0 3.15];
% baseWindow=[2.25 3.5];
% ds=[];
% ds=2;
% ds=2;
% ds=20;
ds=1;
doMed=0;
doFill=1;

if length(trialDuration)>1
    t=trialDuration;
else
%     t=linspace(0,trialDuration,size(data1,2));
    t=linspace(1,trialDuration-1,size(data1,2));
end

% if baseSubtract==1
%     for i=1:size(data1,1)
%         base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
%         data1(i,:)=data1(i,:)-base1;
%         data2(i,:)=data2(i,:)-base2;
%     end
% end

if norm==1
    for i=1:size(data1,1)
        scale=nanmean(data1(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data1(i,:)=data1(i,:)./scale;
%         scale=nanmean(data2(i,t>=normWindow(1) & t<=normWindow(2)),2);
        data2(i,:)=data2(i,:)./scale;
    end
end
  
if baseSubtract==1
    for i=1:size(data1,1)
        base1=nanmean(data1(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        base2=nanmean(data2(i,t>=baseWindow(1) & t<=baseWindow(2)),2);
        data1(i,:)=data1(i,:)-base1;
        data2(i,:)=data2(i,:)-base2;
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

if doMed==1
    if length(trialDuration)>1
        figure();
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
        hold on;
        if doFill==1
            fill([trialDuration fliplr(trialDuration)],[prctile(data1,55) fliplr(prctile(data1,45))],[0.5 0.5 0.5]);
        end
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
%         plot(trialDuration,nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
%         plot(trialDuration,nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        
        plot(trialDuration,nanmedian(data2,1),'Color',c2);
        if doFill==1
            fill([trialDuration fliplr(trialDuration)],[prctile(data2,55) fliplr(prctile(data2,45))],[0.1 0.7 0.5]);
        end
        plot(trialDuration,nanmedian(data2,1),'Color',c2);
%         plot(trialDuration,nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
%         plot(trialDuration,nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        
        plot(trialDuration,nanmedian(data1,1),'Color',c1);
    else
        figure();
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
        hold on;
        if doFill==1
            fill([linspace(0,trialDuration,size(data1,2)) fliplr(linspace(0,trialDuration,size(data1,2)))],[nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)) fliplr(nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)))],[0.5 0.5 0.5]);
        end
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1)+nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1)-nanstd(data1,[],1)./sqrt(size(data1,1)),'Color',c1);
        
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1),'Color',c2);
        if doFill==1
            fill([linspace(0,trialDuration,size(data2,2)) fliplr(linspace(0,trialDuration,size(data2,2)))],[nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)) fliplr(nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)))],[0.1 0.7 0.5]);
        end
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1),'Color',c2);
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1)+nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        plot(linspace(0,trialDuration,size(data2,2)),nanmedian(data2,1)-nanstd(data2,[],1)./sqrt(size(data2,1)),'Color',c2);
        
        plot(linspace(0,trialDuration,size(data1,2)),nanmedian(data1,1),'Color',c1);
    end
    return
end

if length(trialDuration)>1
    dropthebase=[0 4];
    average_baseline_data1=nanmean(nanmean(data1(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
    average_baseline_data2=nanmean(nanmean(data2(:,trialDuration>=dropthebase(1) & trialDuration<=dropthebase(2)),2),1);
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
