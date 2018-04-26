function [withinCondition,acrossCondition,pvalMatrix,pvalForPopvecs,returnPopVec_library,returnPopVec_compare]=comparePopulationVectors(spikes,useTheseAssigns,stimCond_library,ledCond_library,stimCond_compare,ledCond_compare,timeWindow,trials,nBins,popVecL,popVecC)

useAbsoluteValueOfDistances=1;
showFigs=0;
normalizeByVariance=1;

if isempty(useTheseAssigns)
else
    spikes=filtspikes(spikes,0,'assigns',useTheseAssigns);
end

a=unique(spikes.assigns);

% Construct baseline population vector library
% Rows are different units, columns are different trials
if isempty(popVecL)
    library_spikes=filtspikes(spikes,0,'stimcond',stimCond_library,'led',ledCond_library);
    if isempty(trials)
        t=unique(library_spikes.trials);
    else
        t=trials;
    end
    library_trials=t;
    popVec_library=zeros(length(a),length(t));
    for i=1:length(t)
        % Get spike rate of each unit for each trial during timeWindow
        currSpikes_byTrial=filtspikes(library_spikes,0,'trials',t(i));
        for j=1:length(a)
            currSpikes=filtspikes(currSpikes_byTrial,0,'assigns',a(j)); % spikes for current trial for current unit
            numSpikes=sum(currSpikes.spiketimes>timeWindow(1) & currSpikes.spiketimes<timeWindow(2));
            spikeRate=numSpikes/(timeWindow(2)-timeWindow(1));
            popVec_library(j,i)=spikeRate;
        end
    end
else
    popVec_library=popVecL;
end

% Bin pop vec firing rates
% Make nBins for each unit (row)
% binMatrix=zeros(length(a),nBins+1);
% binnedPopVec_library=zeros(length(a),length(t));
% for i=1:size(popVec_library,1)
%     minS=min(popVec_library(i,:));
%     maxS=max(popVec_library(i,:));
%     if minS==maxS
%         % no need to bin -- everything is the same
%     else
%         binSize=(maxS-minS)/nBins;
%         bins=minS:binSize:maxS;
%         binMatrix(i,:)=bins;
%         for j=1:size(popVec_library,2)
%             currSpikeRate=popVec_library(i,j);
%             for k=1:length(bins)-1
%                 if currSpikeRate>=bins(k) && currSpikeRate<=bins(k+1)
%                     binnedPopVec_library(i,j)=(bins(k)+bins(k+1))/2;
%                     break
%                 end
%             end
%         end
%     end
% end
% % Construct bin pdf
% binCounts=zeros(length(a),nBins);
% for i=1:size(binnedPopVec_library,1)
%     for j=1:size(binMatrix,2)-1
%         binCounts(i,j)=sum(binnedPopVec_library(i,:)>=binMatrix(i,j) & binnedPopVec_library(i,:)<binMatrix(i,j+1));
%     end
% end

% Construct comparison population vector library
if isempty(popVecC)
    compare_spikes=filtspikes(spikes,0,'stimcond',stimCond_compare,'led',ledCond_compare);
    if isempty(trials)
        t=unique(compare_spikes.trials);
    else
        t=trials;
    end
    compare_trials=t;
    popVec_compare=zeros(length(a),length(t));
    for i=1:length(t)
        % Get spike rate of each unit for each trial during timeWindow
        currSpikes_byTrial=filtspikes(compare_spikes,0,'trials',t(i));
        for j=1:length(a)
            currSpikes=filtspikes(currSpikes_byTrial,0,'assigns',a(j)); % spikes for current trial for current unit
            numSpikes=sum(currSpikes.spiketimes>timeWindow(1) & currSpikes.spiketimes<timeWindow(2));
            spikeRate=numSpikes/(timeWindow(2)-timeWindow(1));
            popVec_compare(j,i)=spikeRate;
        end
    end
else
    popVec_compare=popVecC;
end
    
% Find visually responsive units - calculate p-value for each unit
% Assume trials are independent
% pvalMatrix=zeros(size(popVec_compare,1),size(popVec_compare,2));
% for i=1:size(popVec_compare,1)
%     for j=1:size(popVec_compare,2)
%         foundMatch=0;
%         for k=1:size(binCounts,2)
%             if popVec_compare(i,j)>=binMatrix(i,k) & popVec_compare(i,j)<=binMatrix(i,k+1)
%                 pval=binCounts(i,k)/sum(binCounts(i,:));
%                 pvalMatrix(i,j)=pval;
%                 foundMatch=1;
%                 break
%             end
%         end
%         if foundMatch==0 || pval==0
%             pvalMatrix(i,j)=1/(sum(binCounts(i,:))+1);
%         end
%     end
% end
pvalMatrix=zeros(size(popVec_compare,1),1);
for i=1:size(popVec_compare,1)
    pvalMatrix(i)=ranksum(popVec_library(i,:),popVec_compare(i,:));
end
        
returnPopVec_library=popVec_library;
returnPopVec_compare=popVec_compare;

% Calculate variance of each row
vars=var([popVec_library popVec_compare]');
% If variance is zero, unit is not informative so throw out this row
popVec_library=popVec_library(vars~=0,:);
popVec_compare=popVec_compare(vars~=0,:);
vars=vars(vars~=0);

useAbsoluteValueDistances=0;
% Calculate pop vec distances within a condition (within library condition)
% Normalize distances by variance of each row
library_distances=[];
for i=1:size(popVec_library,2)
%     for j=i+1:size(popVec_library,2)
    for j=1:size(popVec_library,2)
        if useAbsoluteValueDistances==1
            if normalizeByVariance==1
                d=sum(abs(popVec_library(:,j)-popVec_library(:,i))./vars');
            else
                d=sum(abs(popVec_library(:,j)-popVec_library(:,i)));
            end
        else
            if normalizeByVariance==1
                d=sum((popVec_library(:,j)-popVec_library(:,i))./vars');
            else
                d=sum((popVec_library(:,j)-popVec_library(:,i)));
            end
        end
        library_distances=[library_distances; d];
    end
end

% Calculate pop vec distances across conditions (vs compare condition)
% Normalize distances by variance of each row in library
compare_distances=[];
for i=1:size(popVec_library,2)
    for j=1:size(popVec_compare,2)
        if useAbsoluteValueDistances==1
            if normalizeByVariance==1
                d=sum(abs(popVec_compare(:,j)-popVec_library(:,i))./vars');
            else
                d=sum(abs(popVec_compare(:,j)-popVec_library(:,i)));
            end
        else
            if normalizeByVariance==1
                d=sum((popVec_compare(:,j)-popVec_library(:,i))./vars');
            else
                d=sum((popVec_compare(:,j)-popVec_library(:,i)));
            end
        end
        compare_distances=[compare_distances; d];
    end
end        

library_distances=library_distances(~isnan(library_distances));
compare_distances=compare_distances(~isnan(compare_distances));
if useAbsoluteValueOfDistances==1
    library_distances=abs(library_distances);
    compare_distances=abs(compare_distances);
end
% Plot distributions for pop vec distances within a condition vs compare condition 
pvalForPopvecs=ranksum(library_distances,compare_distances);
[n1,x1]=hist(library_distances,30);
[n2,x2]=hist(compare_distances,30);
if showFigs==1
    figure(); plot(x1,n1,'Color','black'); hold on; plot(x2,n2,'Color','red');
end
% sn1=sum(n1);
% sn2=sum(n2);
% x=sn2/sn1;
[n1,x1]=histnorm(library_distances,30);
[n2,x2]=histnorm(compare_distances,30);
if showFigs==1
    figure(); plot(x1,n1,'Color','black'); hold on; plot(x2,n2,'Color','red');
end
withinCondition.x=x1;
withinCondition.n=n1;
withinCondition.distances=library_distances;
acrossCondition.x=x2;
acrossCondition.n=n2;
acrossCondition.distances=compare_distances;