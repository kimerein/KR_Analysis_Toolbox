function LFPsweepCorrelations_withDistance(LFPbySweep,stimsForSweeps,timeVector)

LFPbySweep=LFPbySweep(~isnan(stimsForSweeps),:);
stimsForSweeps=stimsForSweeps(~isnan(stimsForSweeps));

diffStimConds=sort(unique(stimsForSweeps));
%diffStimConds=diffStimConds(~isnan(diffStimConds));

averages=zeros(length(diffStimConds),size(LFPbySweep,2));
for i=1:length(diffStimConds)
    theseTrials=find(stimsForSweeps==diffStimConds(i));
    takeInds=randperm(length(theseTrials));
    trialsForConds=theseTrials(sort(takeInds(1:67)));
    averages(i,:)=mean(LFPbySweep(trialsForConds,:),1);
end
for i=1:length(diffStimConds)
    averages(i,:)=mean(LFPbySweep(stimsForSweeps==diffStimConds(i),:),1);
end

% Align averages to means
for i=1:length(diffStimConds)
    averages(i,:)=averages(i,:)-mean(averages(i,:));
end

avPointByPointDifference=zeros(length(diffStimConds),length(diffStimConds));
for i=1:length(diffStimConds)
    for j=1:length(diffStimConds)
        avPointByPointDifference(i,j)=mean(abs(averages(i,:)-averages(j,:)));
    end
end

avPointByPointDifference

similarityScore=1-avPointByPointDifference;
similarityScore

% Align trials to means
for i=1:size(LFPbySweep,1)
    LFPbySweep(i,:)=LFPbySweep(i,:)-mean(LFPbySweep(i,:));
end

differences=zeros(size(LFPbySweep,1),size(LFPbySweep,1));
for i=1:size(LFPbySweep,1)
    for j=i:size(LFPbySweep,1)
        pbpDiff=mean(abs(LFPbySweep(i,:)-LFPbySweep(j,:)));
        differences(i,j)=pbpDiff;
    end
end

hold off
figure;
[n,xout]=hist(differences(:),30);
bar(xout,n);

all_ds=zeros(length(diffStimConds),length(diffStimConds));
stdev_all_ds=zeros(length(diffStimConds),length(diffStimConds));
figure;
for i=1:length(diffStimConds)
    for j=i:length(diffStimConds)
        stimCond1=diffStimConds(i);
        stimCond2=diffStimConds(j);
        stim1Inds=find(stimsForSweeps==stimCond1);
        stim2Inds=find(stimsForSweeps==stimCond2);
        d_s=zeros(length(stim1Inds)*length(stim2Inds),1);
        ind=1;
        for k=1:length(stim1Inds)
            for l=1:length(stim2Inds)
                d_s(ind)=differences(stim1Inds(k),stim2Inds(l));
                ind=ind+1;
            end
        end
        all_ds(i,j)=mean(d_s);
        stdev_all_ds(i,j)=std(d_s);
        subplot(length(diffStimConds),length(diffStimConds),(i-1)*length(diffStimConds)+j);
        [xo,no] = histnorm(1-d_s, 30);
        bar(no,xo,'c');
        axis([0.6 1 0 6]);
        %axis 'auto x'        
    end
end

getWithinCond_Ds=zeros(length(diffStimConds),1);
for i=1:length(diffStimConds)
    getWithinCond_Ds(i)=all_ds(i,i);
end
getBetweenCond_Ds=[];
for i=1:length(diffStimConds)
    for j=i+1:length(diffStimConds)
        getBetweenCond_Ds=[getBetweenCond_Ds; all_ds(i,j)];
    end
end
  
ss_withinCond=1-getWithinCond_Ds;
ss_betweenCond=1-getBetweenCond_Ds;

figure;
title('Distributions of Similarity Score (SS) for Within vs. Between Stimulus Conditions - Red (within)');
[xo,no] = histnorm(ss_withinCond, 5);
bar(no,xo,'r');
hold on
[xo,no] = histnorm(ss_betweenCond, 5);
bar(no,xo,'b');


%'r', 'g', 'b', 'c', 'm', 'y', 'k'
        
        
                
                
                