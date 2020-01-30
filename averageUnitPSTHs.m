function [blackPSTH,redPSTH]=averageUnitPSTHs(varargin)

if length(varargin)==1
    psth=varargin{1};
    useStims=[];
elseif length(varargin)==2
    psth=varargin{1};
    useStims=varargin{2};
end

downSamp=1; % down sample factor
% downSamp=10; % down sample factor
blackLED=[0];
% blackLED=[0 0.05];
% freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% blackLED=freqs;
% redLED=[0.05 0.7575 0.2525 0.25 2.5 5.00 5.05 6.06 12.12 50.5];
% redLED=[single(50.5) double(50.5)];
redLED=[0];
% redLED=[double(6.06)];
suppressFigures=0;
% redLED=[0];
% blackLED=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% redLED=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60]+0.05;

% l=psth.unitLED{1};

blackPSTH=zeros(length(psth.psths),length(psth.t));
redPSTH=zeros(length(psth.psths),length(psth.t));
for i=1:length(psth.psths)
    p=psth.psths{i};
    l=psth.unitLED{i};
    s=psth.unitStimcond{i};
    if ~isempty(useStims)
        blackPSTH(i,:)=nanmean(p(ismember(l,blackLED) & ismember(s,useStims),:),1);
        redPSTH(i,:)=nanmean(p(ismember(l,redLED) & ismember(s,useStims),:),1);
    else
        blackPSTH(i,:)=nanmean(p(ismember(l,blackLED),:),1);
        redPSTH(i,:)=nanmean(p(ismember(l,redLED),:),1);
    end
end

if suppressFigures==1
    return
end

% figure(); 
% hax=axes();
downt=downSampAv(psth.t,downSamp);
% downy=downSampMatrix(nanmean(blackPSTH,1),downSamp);
% h=plot(downt,downy,'Color','k');
sub_plotWStderr(downSampMatrix(blackPSTH,downSamp),downSampMatrix(redPSTH,downSamp),downt,'k','r');
% hold on;
% % addErrBar(downt,downy,nanstd(downSampMatrix(blackPSTH,downSamp),[],1)./sqrt(size(blackPSTH,1)),'y',hax,h);
% downy=downSampMatrix(nanmean(redPSTH,1),downSamp);
% h=plot(downt,downy,'Color','r');
% % addErrBar(downt,downy,nanstd(downSampMatrix(redPSTH,downSamp),[],1)./sqrt(size(redPSTH,1)),'y',hax,h);

% blackDist=blackPSTH(:,downt>=3 & downt<=4);
% redDist=redPSTH(:,downt>=3 & downt<=4);
% tDist=downt(downt>=3 & downt<=4);
% for i=1:20
%     currBlack=blackDist(:,i);
%     currRed=redDist(:,i);
%     [h,prank]=ttest2(currBlack,currRed);
% %     disp([prank tDist(i)]);  
% end

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
