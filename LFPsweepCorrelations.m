function LFPsweepCorrelations(LFPbySweep,stimsForSweeps)

LFPbySweep=LFPbySweep(~isnan(stimsForSweeps),:);
stimsForSweeps=stimsForSweeps(~isnan(stimsForSweeps));

diffStimConds=sort(unique(stimsForSweeps));
%diffStimConds=diffStimConds(~isnan(diffStimConds));

averages=zeros(length(diffStimConds),size(LFPbySweep,2));
for i=1:length(diffStimConds)
    averages(i,:)=mean(LFPbySweep(stimsForSweeps==diffStimConds(i),:),1);
end

[av_corrcoefs,av_p,av_rlo,av_rup]=corrcoef(averages');

av_corrcoefs
av_p
av_rlo
av_rup

[corrcoefs,p,rlo,rup]=corrcoef(LFPbySweep');

hold off
figure;
[n,xout]=hist(corrcoefs(:),30);
bar(xout,n);

all_ps=zeros(length(diffStimConds),length(diffStimConds));
all_rs=zeros(length(diffStimConds),length(diffStimConds));
stdev_all_ps=zeros(length(diffStimConds),length(diffStimConds));
stdev_all_rs=zeros(length(diffStimConds),length(diffStimConds));
figure;
for i=1:length(diffStimConds)
    for j=i:length(diffStimConds)
        stimCond1=diffStimConds(i);
        stimCond2=diffStimConds(j);
        stim1Inds=find(stimsForSweeps==stimCond1);
        stim2Inds=find(stimsForSweeps==stimCond2);
%         p_sum=0;
%         r_sum=0;
%         p_stdev=0;
%         r_stdev=0;
        p_s=zeros(length(stim1Inds)*length(stim2Inds),1);
        r_s=zeros(length(stim1Inds)*length(stim2Inds),1);
        ind=1;
        for k=1:length(stim1Inds)
            for l=1:length(stim2Inds)
%                 r_sum=r_sum+corrcoefs(stim1Inds(k),stim2Inds(l));
%                 p_sum=p_sum+p(stim1Inds(k),stim2Inds(l));
                p_s(ind)=p(stim1Inds(k),stim2Inds(l));
                r_s(ind)=corrcoefs(stim1Inds(k),stim2Inds(l));
                ind=ind+1;
            end
        end
        all_ps(i,j)=mean(p_s);
        all_rs(i,j)=mean(r_s);
        stdev_all_ps(i,j)=std(p_s);
        stdev_all_rs(i,j)=std(r_s);
        subplot(length(diffStimConds),length(diffStimConds),(i-1)*length(diffStimConds)+j);
        [xo,no] = histnorm(r_s, 30);
        bar(no,xo,'c');
    end
end

getWithinCond_Ps=zeros(length(diffStimConds),1);
for i=1:length(diffStimConds)
    getWithinCond_Ps(i)=all_ps(i,i);
end
getBetweenCond_Ps=[];
for i=1:length(diffStimConds)
    for j=i+1:length(diffStimConds)
        getBetweenCond_Ps=[getBetweenCond_Ps; all_ps(i,j)];
    end
end
getWithinCond_Rs=zeros(length(diffStimConds),1);
for i=1:length(diffStimConds)
    getWithinCond_Rs(i)=all_rs(i,i);
end
getBetweenCond_Rs=[];
for i=1:length(diffStimConds)
    for j=i+1:length(diffStimConds)
        getBetweenCond_Rs=[getBetweenCond_Rs; all_rs(i,j)];
    end
end
         
figure;
title('Distributions of R for Within vs. Between Stimulus Conditions - Red (within)');
[xo,no] = histnorm(getWithinCond_Rs, 30);
bar(no,xo,'r');
hold on
[xo,no] = histnorm(getBetweenCond_Rs, 30);
bar(no,xo,'b');

figure;
title('Distributions of P for Within vs. Between Stimulus Conditions - Red (within)');
[xo,no] = histnorm(getWithinCond_Ps, 30);
bar(no,xo,'r');
hold on
[xo,no] = histnorm(getBetweenCond_Ps, 30);
bar(no,xo,'b');


%'r', 'g', 'b', 'c', 'm', 'y', 'k'
        
        
                
                
                