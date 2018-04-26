function [temp_c,allcin1,allcin2,allcinRef]=analyzeTransitionCorrelation(currs,backup_allSfirst,backup_allSsecond,backup_ffr,finterp,currBand,ffrLEDadjusted2,evBase)

sumcells=0;

allSfirst=backup_allSfirst;
allSsecond=backup_allSsecond;
ffrLEDadjusted=backup_ffr;

if ~iscell(ffrLEDadjusted)
    % Subtract means off all units
    for i=1:length(allSfirst)
        %     currUnitMin=nanmin(nanmin([allSfirst{i} allSsecond{i}],[],2),[],1);
        %     allSfirst{i}=allSfirst{i}-currUnitMin;
        %     allSsecond{i}=allSsecond{i}-currUnitMin;
        %     currUnitMax=nanmax(nanmax([allSfirst{i} allSsecond{i}],[],2),[],1);
        %     allSfirst{i}=allSfirst{i}./currUnitMax;
        %     allSsecond{i}=allSsecond{i}./currUnitMax;
        %     ffrLEDadjusted(i,:)=ffrLEDadjusted(i,:)-min(ffrLEDadjusted(i,:));
        %     ffrLEDadjusted(i,:)=ffrLEDadjusted(i,:)./max(ffrLEDadjusted(i,:));
        currUnitMean=nanmean(nanmean([allSfirst{i} allSsecond{i}],2),1);
%         allSfirst{i}=allSfirst{i}-currUnitMean;
%         allSsecond{i}=allSsecond{i}-currUnitMean;
        %     ffrLEDadjusted(i,:)=ffrLEDadjusted(i,:)-currUnitMean;
        allSfirst{i}=allSfirst{i}-repmat(ffrLEDadjusted(i,:)',1,size(allSfirst{i},2));
        allSsecond{i}=allSsecond{i}-repmat(ffrLEDadjusted(i,:)',1,size(allSfirst{i},2));
    end
else
    for i=1:length(allSfirst)
%         allSfirst{i}=allSfirst{i}-(1/4)*ffrLEDadjusted{i};
%         allSsecond{i}=allSsecond{i}-(1/4)*ffrLEDadjusted{i};

%         % Using
%         allSfirst{i}=allSfirst{i}-ffrLEDadjusted{i};
%         allSsecond{i}=allSsecond{i}-ffrLEDadjusted{i};
        
%         allSfirst{i}=allSfirst{i}-spontBase(i);
%         allSsecond{i}=allSsecond{i}-spontBase(i);

%         % Using
%         currUnitMean=nanmean(nanmean([allSfirst{i} allSsecond{i}],2),1);
%         allSfirst{i}=allSfirst{i}-currUnitMean;
%         allSsecond{i}=allSsecond{i}-currUnitMean;

%         % Using new
%         currUnitMean=nanmean([allSfirst{i}; allSsecond{i}],1);
%         allSfirst{i}=allSfirst{i}-repmat(currUnitMean,12,1);
%         allSsecond{i}=allSsecond{i}-repmat(currUnitMean,12,1);
        
%         temp=allSfirst{i};
%         for j=1:size(temp,1)
%             temp(j,:)=temp(j,:)-min(temp(j,:));
%             temp(j,:)=temp(j,:)./max(temp(j,:));
%         end
%         allSfirst{i}=temp;
%         temp=allSsecond{i};
%         for j=1:size(temp,1)
%             temp(j,:)=temp(j,:)-min(temp(j,:));
%             temp(j,:)=temp(j,:)./max(temp(j,:));
%         end
%         allSsecond{i}=temp;
        
%         allSfirst{i}=allSfirst{i}+evBase(i)./14.5; % Convert from spikes to freq amp
%         allSsecond{i}=allSsecond{i}+evBase(i)./14.5;
        
%         % Using
%         allSfirst{i}=allSfirst{i}+evBase(i)./60; % Convert from spikes to freq amp
%         allSsecond{i}=allSsecond{i}+evBase(i)./60;
    end
end
        
if sumcells==1
    if currs==12
        prefSet1=find(ffrLEDadjusted2(:,1)>=prctile(ffrLEDadjusted2(:,1),75));
        prefSet2=find(ffrLEDadjusted2(:,1)>=prctile(ffrLEDadjusted2(:,1),50) & ffrLEDadjusted2(:,1)<prctile(ffrLEDadjusted2(:,1),75));
        nonPrefSet1=find(ffrLEDadjusted2(:,1)<=prctile(ffrLEDadjusted2(:,1),25));
        nonPrefSet2=find(ffrLEDadjusted2(:,1)>prctile(ffrLEDadjusted2(:,1),25) & ffrLEDadjusted2(:,1)<prctile(ffrLEDadjusted2(:,1),50)); 
    else
        prefSet1=find(ffrLEDadjusted2(:,currs+1)>=prctile(ffrLEDadjusted2(:,currs+1),75));
        prefSet2=find(ffrLEDadjusted2(:,currs+1)>=prctile(ffrLEDadjusted2(:,currs+1),50) & ffrLEDadjusted2(:,currs+1)<prctile(ffrLEDadjusted2(:,currs+1),75));
        nonPrefSet1=find(ffrLEDadjusted2(:,currs+1)<=prctile(ffrLEDadjusted2(:,currs+1),25));
        nonPrefSet2=find(ffrLEDadjusted2(:,currs+1)>prctile(ffrLEDadjusted2(:,currs+1),25) & ffrLEDadjusted2(:,currs+1)<prctile(ffrLEDadjusted2(:,currs+1),50)); 
    end
    temp_prefSet1=zeros(size(allSfirst{1}));
    temp_nonPrefSet1=zeros(size(allSfirst{1}));
    temp_prefSet2=zeros(size(allSfirst{1}));
    temp_nonPrefSet2=zeros(size(allSfirst{1}));
    for i=1:length(allSfirst)
        if ismember(i,prefSet1)
            temp_prefSet1=temp_prefSet1+allSfirst{i};
        elseif ismember(i,prefSet2)
            temp_prefSet2=temp_prefSet2+allSfirst{i};
        elseif ismember(i,nonPrefSet1)
            temp_nonPrefSet1=temp_nonPrefSet1+allSfirst{i};
        elseif ismember(i,nonPrefSet2)
            temp_nonPrefSet2=temp_nonPrefSet2+allSfirst{i};
        end
    end
    clear allSfirst 
    allSfirst{1}=temp_prefSet1;
    allSfirst{2}=temp_prefSet2;
    allSfirst{3}=temp_nonPrefSet1;
    allSfirst{4}=temp_nonPrefSet2;
    temp_prefSet1=zeros(size(allSsecond{1}));
    temp_nonPrefSet1=zeros(size(allSsecond{1}));
    temp_prefSet2=zeros(size(allSsecond{1}));
    temp_nonPrefSet2=zeros(size(allSsecond{1}));
    for i=1:length(allSsecond)
        if ismember(i,prefSet1)
            temp_prefSet1=temp_prefSet1+allSsecond{i};
        elseif ismember(i,prefSet2)
            temp_prefSet2=temp_prefSet2+allSsecond{i};
        elseif ismember(i,nonPrefSet1)
            temp_nonPrefSet1=temp_nonPrefSet1+allSsecond{i};
        elseif ismember(i,nonPrefSet2)
            temp_nonPrefSet2=temp_nonPrefSet2+allSsecond{i};
        end
    end
    clear allSsecond
    allSsecond{1}=temp_prefSet1;
    allSsecond{2}=temp_prefSet2;
    allSsecond{3}=temp_nonPrefSet1;
    allSsecond{4}=temp_nonPrefSet2;
    unitsLike=1:4;
end
            

% Choose units that like this stim
% currs=1;
% unitsLike=find(ffrLEDadjusted(:,currs)>nanmedian(ffrLEDadjusted(:,currs)));
% unitsLike=find(ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50));
% unitsLike=find(evBase(:,currs)<prctile(evBase(:,currs),50));
% unitsLike=find(evBase(:,currs)>0 & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50));
% unitsLike=find(evBase(:,currs)>prctile(evBase(:,currs),30) & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50));
% unitsLike=find(evBase(:,currs)>prctile(evBase(:,currs),40) & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50));
% unitsLike=find(~(evBase(:,currs)<prctile(evBase(:,currs),70) & ffrLEDadjusted2(:,currs)<prctile(ffrLEDadjusted2(:,currs),70)));
% unitsLike=find((evBase(:,currs)>prctile(evBase(:,currs),40) & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),30)) & ~(evBase(:,currs)<prctile(evBase(:,currs),60) & ffrLEDadjusted2(:,currs)<prctile(ffrLEDadjusted2(:,currs),50)));
% unitsLike=find(~(evBase(:,currs)<prctile(evBase(:,currs),50) & ffrLEDadjusted2(:,currs)<prctile(ffrLEDadjusted2(:,currs),50)));
% unitsLike=find(~(evBase(:,currs)<prctile(evBase(:,currs),70) & ffrLEDadjusted2(:,currs)<prctile(ffrLEDadjusted2(:,currs),50)));
% unitsLike=find(~(evBase(:,currs)>prctile(evBase(:,currs),40) & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50)));
% unitsLike=find(evBase(:,currs)>prctile(evBase(:,currs),40) & ffrLEDadjusted2(:,currs)>prctile(ffrLEDadjusted2(:,currs),50));
e=evBase(:,currs);
[~,evInd]=sort(e);
[~,ffrInd]=sort(ffrLEDadjusted2(:,currs));
combInd=zeros(1,length(evInd));
for i=1:length(evInd)
    combInd(i)=nanmean([find(i==evInd) find(i==ffrInd)]);
end
getn=4;
[~,sortInd]=sort(combInd);
unitsLike=sortInd(end-(getn-1):end);
unitsLike=evInd(end-(getn-1):end);
% unitsLike=ffrInd(end-(getn-1):end);
unitsLike=1:13;

% [~,evInd]=sort(-e);
% unitsOppLike=evInd(end-(getn-1):end);
% [unitsLike,oppInd]=sort([unitsLike; unitsOppLike]);
% areOpp=ismember(oppInd,getn+1:getn+getn);
areOpp=zeros(1,length(unitsLike));

% if currs==12
%     unitsLike=find(evBase(:,1)>prctile(evBase(:,1),40) & ffrLEDadjusted2(:,1)>prctile(ffrLEDadjusted2(:,1),50));
% else
%     unitsLike=find(evBase(:,currs+1)>prctile(evBase(:,currs+1),40) & ffrLEDadjusted2(:,currs+1)>prctile(ffrLEDadjusted2(:,currs+1),50));
% end

% unitsLike=find(ffrLEDadjusted2(:,currs)>nanmedian(ffrLEDadjusted2(:,currs)));
% if currs==12
%     unitsLike=find(ffrLEDadjusted(:,1)<nanmedian(ffrLEDadjusted(:,1)));
% else
%     unitsLike=find(ffrLEDadjusted(:,currs)<nanmedian(ffrLEDadjusted(:,currs)));
% end
% unitsLike=1:13;

% Analyze transition correlation
freqBands={[3 20];[3 3.5];[3.5 5];[5.5 6.5];[6 7.5];[7.5 8.5];[8.5 9.5];[9.5 10.5];[10.5 11.5];[11.5 12.5];[12.5 13.5];[13.5 14.5];[14.5 20]};
for i=1:length(unitsLike)
temp=allSfirst{unitsLike(i)};
% if currs==12
%     currStimAndFreq(i)=evBase(i,1);
% else
%     currStimAndFreq(i)=evBase(i,currs+1);
% end
% currStimAndFreq(i)=evBase(i,currs);
currStimAndFreq(i)=nanmean(temp(currs,finterp>=currBand(1) & finterp<=currBand(2)),2);
temp=allSsecond{unitsLike(i)}; % for persistence
% currStimAndFreqSameTrial(i)=nanmean(temp(currs,finterp>=currBand(1) & finterp<=currBand(2)),2);
currStimAndFreqSameTrial(i)=nanmean(temp(currs,finterp>=8.5 & finterp<=9.5),2);
% if currs==12
%     currStimAndFreqSameTrial(i)=evBase(i,1);
% else
%     currStimAndFreqSameTrial(i)=evBase(i,currs+1);
% end
end
[~,si]=sort(currStimAndFreq);
cin1=currStimAndFreq(si);
newOpp=areOpp(si);
oi=find(newOpp==1);
% cin1=-cin1;
cinRef=currStimAndFreqSameTrial(si);
% cinRef=-cinRef;
allcin1=cin1;
allcinRef=cinRef;
% figure(); % PLOT
% subplot(length(freqBands)+1,1,1); % PLOT
% % heatmap(currStimAndFreq(si)); % PLOT
% heatmap(cin1); % PLOT
allcin2=[];
for j=1:length(freqBands)
currFreqBand=freqBands{j};
for i=1:length(unitsLike)
temp=allSsecond{unitsLike(i)};
currStimAndFreq2(i)=nanmean(temp(currs,finterp>=currFreqBand(1) & finterp<=currFreqBand(2)),2);
end
% subplot(length(freqBands)+1,1,j+1); % PLOT
% heatmap(currStimAndFreq2(si)); % PLOT
cin2=currStimAndFreq2(si);
cin2(oi)=cin2(fliplr(oi'));
% if all(cin1==cinRef)
%     c=nan;
%     p=0;
% else
%     [b,bint,r,rint,stats]=regress(cin2',[ones(size(cin1')) cin1' cinRef']); % for persistence
%     % cin2=cin2-b(3).*cinRef;
%     if any(b<0)
%         c=0;
%     else
%         c=b(2)/(b(2)+b(3));
%     end
%     % % c=b(2);
%     p=stats(3);
%     if p>0.05
%         c=nan;
%     end
% end
% if p>0.05
%     c=nan;
% end
% temp_c(j)=c;

% [c,p]=corr(cin1',cin2','type','Spearman');
[ctemp,ptemp]=corrcoef(cin1,cin2);
% [cref,pref]=corrcoef(cinRef,cin2);
% c=ctemp(1,2)-cref(1,2);
p=ptemp(1,2);
c=ctemp(1,2);
% if ptemp>0.3
if ptemp>1
    c=nan;
end
% if p>0.05
%     c=0;
% end
% c=cref(1,2);
temp_c(j)=c;
allcin2=[allcin2; cin2];
% disp([freqBands{j} c p]); % PLOT 
end
% figure(); 
% scatter(allcin1,allcin2(11,:));
% figure(); 
% subplot(2,1,1);
% heatmap(allcin1);
% subplot(2,1,2);
% heatmap(allcin2(11,:));
% realtrans_c(1,:)=temp_c;