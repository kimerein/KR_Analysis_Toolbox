function [stim1rawsum,stim2rawsum,stim3rawsum,psths_trialByTrial,prefStimset,psths_ref,sigStim,arrayVecs]=doPopVecPredictionAnalysis(grp1_spikes,grp2_spikes,allAssigns,alphaYes,gp1_unitByUnit,gp2_unitByUnit,useStimcond,useLED,allRefSpikes)

grp1_responseWindow=[0.3 1.8];
% grp1_responseWindow=[0 1.5];
bin=1; % in ms, 25 is good
trialDuration=3.5; % in s
getTrialByTrial=1;

% Grp 1 -- get unit weight vectors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes=grp1_spikes;

% Get stim. conditions defining different stim. types
if isempty(useStimcond)
    useStimcond{1}=1:4;
    useStimcond{2}=5:8;
    useStimcond{3}=9:12;
end

% Get trials for each unit
if ~isfield(gp1_unitByUnit,'unitByUnitTrials')
    for i=1:length(allAssigns)
        unitByUnitTrials{i}=unique(spikes.sweeps.trials);
        unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
        unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
    end
else
    unitByUnitTrials=gp1_unitByUnit.unitByUnitTrials;
    unitByUnitStimcond=gp1_unitByUnit.unitByUnitStimcond;
    unitByUnitLED=gp1_unitByUnit.unitByUnitLED;
end

% Use alpha spikes or no alpha spikes
% if ~isempty(alphaYes)
%     spikes=filtspikes(spikes,0,'alpha',alphaYes);
% end

% Get unit weight vectors
[unitVecs,prefStimset,sigStim]=getUnitWeightVector(spikes,allAssigns,useStimcond,grp1_responseWindow,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,1,1,useLED);
arrayVecs=nan(length(unitVecs),length(unitVecs{1}));
for i=1:length(unitVecs)
    arrayVecs(i,:)=unitVecs{i};
end

% Grp 2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes=grp2_spikes;

% Get stim. conditions defining different stim. types
if isempty(useStimcond)
    useStimcond{1}=[1 4 7 10];
    useStimcond{2}=[2 5 8 11];
    useStimcond{3}=[3 6 9 12];
end

% Get trials for each unit
if ~isfield(gp2_unitByUnit,'unitByUnitTrials')
    for i=1:length(allAssigns)
        unitByUnitTrials{i}=unique(spikes.sweeps.trials);
        unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
        unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
    end
else
    unitByUnitTrials=gp2_unitByUnit.unitByUnitTrials;
    unitByUnitStimcond=gp2_unitByUnit.unitByUnitStimcond;
    unitByUnitLED=gp2_unitByUnit.unitByUnitLED;
end

% Use alpha spikes or no alpha spikes
if ~isempty(alphaYes)
    spikes=filtspikes(spikes,0,'alpha',alphaYes);
end

% Get PSTH for units
currTrialCon=unitByUnitTrials{1};
currStimCon=unitByUnitStimcond{1};
currLEDCon=unitByUnitLED{1};
if ~isempty(useLED)
    unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{1}) & ismember(currLEDCon,useLED));
else
    unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{1}));
end
[~,~,~,x,y]=psth_wStd(filtspikes(spikes,0,'stimcond',useStimcond{1}),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
psths=nan(length(allAssigns),length(useStimcond),length(x));
psths_t=cell(length(allAssigns),length(useStimcond));
psths_r=cell(1,length(useStimcond));
xt=[];
xtref=[];
for i=1:length(allAssigns)
    for j=1:length(useStimcond)
        useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        [~,~,~,x,psths(i,j,:)]=psth_wStd(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        if getTrialByTrial==1
            [~,~,~,xt,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        end
    end
end
if ~isempty(allRefSpikes)
    disp('ref spikes');
    for j=1:length(useStimcond)  
        disp(j);
        useSpikes=filtspikes(spikes,0,'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{1};
        currStimCon=unitByUnitStimcond{1};
        currLEDCon=unitByUnitLED{1};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        [~,~,~,xtref,~,psths_r{j}]=psth_wStd_trialByTrial(allRefSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
    end
end
psths_trialByTrial.xt=xt;
psths_trialByTrial.psths_t=psths_t;
psths_ref.xt=xtref;
psths_ref.psths_t=psths_r;

% Re-organize psths by current stim.
togethery=nan(length(allAssigns),length(x(x>=0 & x<=1.5))*3);
x1=x(1):x(2)-x(1):4.5;
y=nan(3,length(x(x>=0 & x<=1.5))*3);
for i=1:length(allAssigns)
    temp=makeVecsSameLength(y(1,x1>=0 & x1<=3.5),psths(i,1,:));
    y(1,x1>=0 & x1<=3.5)=temp;
    temp=makeVecsSameLength(y(2,x1>=0 & x1<=0.5),psths(i,2,x>=3));
    y(2,x1>=0 & x1<=0.5)=temp;
    temp=makeVecsSameLength(y(2,find(x1>=1.5,1,'first'):end),psths(i,2,x>=0 & x<=3));
    y(2,find(x1>=1.5,1,'first'):end)=temp;
    temp=makeVecsSameLength(y(3,x1>=0 & x1<=2),psths(i,3,x>=1.5 & x<=3.5));
    y(3,x1>=0 & x1<=2)=temp;
    temp=makeVecsSameLength(y(3,find(x1>=3,1,'first'):end),psths(i,3,x>=0 & x<=1.5));
    y(3,find(x1>=3,1,'first'):end)=temp;
    togethery(i,:)=nanmean(y,1);
end
   
% Add up spikes from unit populations with different preferences
stim1rawsum=nan(length(allAssigns),size(togethery,2));
stim2rawsum=nan(length(allAssigns),size(togethery,2));
stim3rawsum=nan(length(allAssigns),size(togethery,2));  
for i=1:length(allAssigns)
    if sigStim(i)<0.0001
        if prefStimset(i)==1 
%         if rand(1)<0.33
            stim1rawsum(i,:)=togethery(i,:);
%         end
        elseif prefStimset(i)==2
%         if rand(1)<0.33
            stim2rawsum(i,:)=togethery(i,:);
%         end
        elseif prefStimset(i)==3
%         if rand(1)<0.33
            stim3rawsum(i,:)=togethery(i,:);
        end
    end
end
% Norm. each unit PSTH
for i=1:length(allAssigns)
    if ~isnan(stim1rawsum(i,:))
        stim1rawsum(i,:)=stim1rawsum(i,:)-min(stim1rawsum(i,:));
        stim1rawsum(i,:)=stim1rawsum(i,:)./max(stim1rawsum(i,:));
    end
    if ~isnan(stim2rawsum(i,:))
        stim2rawsum(i,:)=stim2rawsum(i,:)-min(stim2rawsum(i,:));
        stim2rawsum(i,:)=stim2rawsum(i,:)./max(stim2rawsum(i,:));
    end
    if ~isnan(stim3rawsum(i,:))
        stim3rawsum(i,:)=stim3rawsum(i,:)-min(stim3rawsum(i,:));
        stim3rawsum(i,:)=stim3rawsum(i,:)./max(stim3rawsum(i,:));
    end
end
normToMean=1;
% Norm. each PSTH section to mean
xtempy{1}=x1>=0 & x1<=1.5;
xtempy{2}=x1>=1.5 & x1<=3;
% xtempy{3}=x1>=3 & x1<=4.5;
xtempy{3}=find(x1>=3,1,'first'):size(stim1rawsum,2);
if normToMean==1
    for i=1:length(xtempy)
        xtemp=xtempy{i};
        stim1rawsum(:,xtemp)=stim1rawsum(:,xtemp)-nanmean(nanmean(stim1rawsum(:,xtemp),2),1);
        stim2rawsum(:,xtemp)=stim2rawsum(:,xtemp)-nanmean(nanmean(stim2rawsum(:,xtemp),2),1);
        stim3rawsum(:,xtemp)=stim3rawsum(:,xtemp)-nanmean(nanmean(stim3rawsum(:,xtemp),2),1);
    end
end
    
figure(); 
hax=axes();
hl=plotLineAndErr(x1,stim1rawsum,'k',hax);
hold on; 
hl=plotLineAndErr(x1,stim2rawsum,'r',hax);
hl=plotLineAndErr(x1,stim3rawsum,'g',hax);
title('Raw Summed Spikes for Units Preferring Each Stim');

% Make representation PSTH
stimVsStim_PSTH=nan(nchoosek(length(useStimcond),2),size(togethery,2));
l=1;
stim1psth=nan(2,size(togethery,2));
stim2psth=nan(2,size(togethery,2));
stim3psth=nan(2,size(togethery,2));
for i=1:length(useStimcond)
    for j=i+1:length(useStimcond)
        stimVsStim_PSTH(l,:)=sum(togethery.*repmat(arrayVecs(:,l),1,size(togethery,2)),1);
        if i==1 && j==2
            a=arrayVecs(:,l);
%             a(a<=0)=0;
            stim1psth(1,:)=sum(togethery.*repmat(a,1,size(togethery,2)),1);
            a=arrayVecs(:,l);
%             a(a>=0)=0;
            stim2psth(1,:)=sum(togethery.*repmat(-a,1,size(togethery,2)),1);
        elseif i==1 && j==3
            a=arrayVecs(:,l);
%             a(a<=0)=0;
            stim1psth(2,:)=sum(togethery.*repmat(a,1,size(togethery,2)),1);
            a=arrayVecs(:,l);
%             a(a>=0)=0;
            stim3psth(1,:)=sum(togethery.*repmat(-a,1,size(togethery,2)),1);
        elseif i==2 && j==3
            a=arrayVecs(:,l);
%             a(a<=0)=0;
            stim2psth(2,:)=sum(togethery.*repmat(a,1,size(togethery,2)),1);
            a=arrayVecs(:,l);
%             a(a>=0)=0;
            stim3psth(2,:)=sum(togethery.*repmat(-a,1,size(togethery,2)),1);
        end
        l=l+1;
    end
end

figure(); 
plot(x1(1:size(stim1psth,2)),nanmean(stim1psth,1),'Color','k');
title('Stim PSTH Weighted by Unit Vec');
hold on; 
plot(x1(1:size(stim1psth,2)),nanmean(stim2psth,1),'Color','r');
plot(x1(1:size(stim1psth,2)),nanmean(stim3psth,1),'Color','g');

figure(); 
plot(x1(1:size(stimVsStim_PSTH,2)),stimVsStim_PSTH(1,:));
title('Stim 1 vs Stim 2');
figure(); 
plot(x1(1:size(stimVsStim_PSTH,2)),stimVsStim_PSTH(2,:));
title('Stim 1 vs Stim 3');
figure(); 
plot(x1(1:size(stimVsStim_PSTH,2)),stimVsStim_PSTH(3,:));
title('Stim 2 vs Stim 3');

% Current vs Previous
cvp(1,:)=stimVsStim_PSTH(2,x1>=0 & x1<=1.5);
te=find(x1>=1.5 & x1<=3);
cvp(2,:)=stimVsStim_PSTH(1,te(1:size(cvp,2)));
te=find(x1>=3,1,'first'):size(stimVsStim_PSTH,2);
cvp(3,:)=stimVsStim_PSTH(3,te);
weighted_curr(1,:)=nanmean(stim1psth(:,x1>=0 & x1<=1.5),1);
weighted_curr(2,:)=nanmean(stim2psth(:,x1>=1.5 & x1<=3),1);
weighted_curr(3,:)=nanmean(stim3psth(:,x1>=3 & x1<=4.5),1);
raw_curr=stim1rawsum(:,x1>=0 & x1<=1.5);
raw_curr=[raw_curr; stim2rawsum(:,x1>=1.5 & x1<=3)];
raw_curr=[raw_curr; stim3rawsum(:,x1>=3 & x1<=4.5)];

% Current vs Next
cvn(1,:)=stimVsStim_PSTH(1,x1>=0 & x1<=1.5);
cvn(2,:)=stimVsStim_PSTH(3,x1>=1.5 & x1<=3);
cvn(3,:)=stimVsStim_PSTH(2,x1>=3 & x1<=4.5);
weighted_next(1,:)=nanmean(stim3psth(:,x1>=0 & x1<=1.5),1);
weighted_next(2,:)=nanmean(stim1psth(:,x1>=1.5 & x1<=3),1);
weighted_next(3,:)=nanmean(stim2psth(:,x1>=3 & x1<=4.5),1);
raw_next=stim3rawsum(:,x1>=0 & x1<=1.5);
raw_next=[raw_next; stim1rawsum(:,x1>=1.5 & x1<=3)];
raw_next=[raw_next; stim2rawsum(:,x1>=3 & x1<=4.5)];

% Previous vs Next
pvn(1,:)=stimVsStim_PSTH(3,x1>=0 & x1<=1.5);
pvn(2,:)=stimVsStim_PSTH(2,x1>=1.5 & x1<=3);
pvn(3,:)=stimVsStim_PSTH(1,x1>=3 & x1<=4.5);
weighted_prev(1,:)=nanmean(stim2psth(:,x1>=0 & x1<=1.5),1);
weighted_prev(2,:)=nanmean(stim3psth(:,x1>=1.5 & x1<=3),1);
weighted_prev(3,:)=nanmean(stim1psth(:,x1>=3 & x1<=4.5),1);
raw_prev=stim2rawsum(:,x1>=0 & x1<=1.5);
raw_prev=[raw_prev; stim3rawsum(:,x1>=1.5 & x1<=3)];
raw_prev=[raw_prev; stim1rawsum(:,x1>=3 & x1<=4.5)];

figure(); 
hax=axes();
hl=plotLineAndErr(x1(x1>=0 & x1<=1.5),raw_curr,'k',hax);
hold on; 
hl=plotLineAndErr(x1(x1>=0 & x1<=1.5),raw_next,'r',hax);
hl=plotLineAndErr(x1(x1>=0 & x1<=1.5),raw_prev,'g',hax);
title('Raw Summed Spikes for Units Preferring Each Stim -- Together');


figure(); 
plot(x1(x1>=0 & x1<=1.5),nanmean(weighted_curr,1),'Color','k');
hold on;
plot(x1(x1>=0 & x1<=1.5),nanmean(weighted_next,1),'Color','r');
plot(x1(x1>=0 & x1<=1.5),nanmean(weighted_prev,1),'Color','g');
legend('curr','next','prev');
title('Weighted');

figure(); 
plot(x1(x1>=0 & x1<=1.5),nanmean(cvp,1),'Color','g');
x_new=x1(x1>=0 & x1<=1.5);
y_new(1,:)=nanmean(cvp,1);
hold on;
plot(x1(x1>=0 & x1<=1.5),nanmean(cvn,1),'Color','r');
y_new(2,:)=nanmean(cvn,1);
plot(x1(x1>=0 & x1<=1.5),nanmean(pvn,1),'Color',[0.5 0.5 0.5]);
y_new(3,:)=nanmean(pvn,1);
legend('curr v prev','curr v next','prev v next');

% Normalize
y_new=y_new./repmat(max(y_new,[],2),1,size(y_new,2));
curr_vs_prev=asin(abs(y_new(1,:)));
curr_vs_prev(y_new(1,:)<0)=-curr_vs_prev(y_new(1,:)<0);
curr_vs_prev=(pi/4)-curr_vs_prev;
curr_vs_next=asin(abs(y_new(2,:)));
curr_vs_next(y_new(2,:)<0)=-curr_vs_next(y_new(2,:)<0);
curr_vs_next=(pi/4)-curr_vs_next;
prev_vs_next=asin(abs(y_new(3,:)));
prev_vs_next(y_new(3,:)<0)=-prev_vs_next(y_new(3,:)<0);
prev_vs_next=(pi/4)-prev_vs_next;

% Estimate population location
est1_x=cos(curr_vs_prev);
est2_x=cos(curr_vs_next);
est1_y=cos(prev_vs_next);
est2_y=sin(curr_vs_prev);
est1_z=sin(prev_vs_next);
est2_z=sin(curr_vs_next);

figure(); 
plot3(nanmean([est1_x(1:floor(end/2)); est2_x(1:floor(end/2))],1),nanmean([est1_y(1:floor(end/2)); est2_y(1:floor(end/2))],1),nanmean([est1_z(1:floor(end/2)); est2_z(1:floor(end/2))],1),'Color','r');
hold on;
plot3(nanmean([est1_x(floor(end/2):end); est2_x(floor(end/2):end)],1),nanmean([est1_y(floor(end/2):end); est2_y(floor(end/2):end)],1),nanmean([est1_z(floor(end/2):end); est2_z(floor(end/2):end)],1),'Color','b');
title('Trajectory in Current (x) Previous (y) Next (z) Space');
xlabel('Current');
ylabel('Previous');
zlabel('Next');
end

function hl=plotLineAndErr(x,data,c,hax)

hl=plot(x(1:size(data,2)),nanmean(data,1),'Color',c);
addErrBar(x(1:size(data,2)),nanmean(data,1),nanstd(data,[],1)./sqrt(size(data,1)),'y',hax,hl);

end

function temp=makeVecsSameLength(a,b)

% Makes vec b same length as vec a

temp=nan(1,length(a));
if length(a)>length(b)
    temp(1:length(b))=b;
elseif length(b)>length(a)
    temp=b(1:length(a));
else
    temp=b;
end

end

function [unitVecs,maxPrefStimset,sigStim]=getUnitWeightVector(spikes,useTheseAssigns,useTheseStimcond,onsetResponseWindow,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,nBest,nWorst,useLED)

prefR=zeros(length(useTheseAssigns),1);
prefStimcond=cell(length(useTheseAssigns),1);
nonprefR=zeros(length(useTheseAssigns),1);
nonprefStimcond=cell(length(useTheseAssigns),1);
sigStim=zeros(length(useTheseAssigns),1);
unitVecs=cell(length(useTheseAssigns),1);
maxPrefStimset=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    stimcondR=zeros(length(useTheseStimcond),1);
    d=cell(1,length(useTheseStimcond));
    for j=1:length(useTheseStimcond)
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useTheseStimcond{j}) & ismember(currLEDCon,useLED));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useTheseStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('unitByUnitConsensus is empty');
        end
        [stimcondR(j),~,d{j}]=calcMeanAndStdForUnit(filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',useTheseStimcond{j}),onsetResponseWindow,unitByUnitConsensus);
    end
    [sortVals,sortInd]=sort(stimcondR,'descend');
    prefR(i)=mean(sortVals(1:nBest));
    ind=sortInd(1:nBest);
    l=[];
    for j=1:length(ind)
        l=[l useTheseStimcond{ind}];
    end
    prefStimcond{i}=l;
%     maxprefStimcond=useTheseStimcond{ind(1)};
    maxprefStimcond=ind(1);
    maxPrefStimset(i)=maxprefStimcond;
    [sortVals,sortInd]=sort(stimcondR,'ascend');
    nonprefR(i)=mean(sortVals(1:nWorst));
    ind=sortInd(1:nWorst);
    l=[];
    for j=1:length(ind)
        l=[l useTheseStimcond{ind}];
    end
    nonprefStimcond{i}=l;
%     minnonprefStimcond=useTheseStimcond{ind(1)};
    minnonprefStimcond=ind(1);
%     sigStim(i)=mattest(d{maxprefStimcond}',d{minnonprefStimcond}');
    sigStim(i)=ranksum(d{maxprefStimcond}',d{minnonprefStimcond}');
    unitWeightVector=zeros(nchoosek(length(useTheseStimcond),2),1);
    l=1;
    for j=1:length(useTheseStimcond)
        for k=j+1:length(useTheseStimcond)
            p=ranksum(d{j},d{k});
            if p>=1
                l=l+1;
                continue
            end
            if nanmean(d{j})>nanmean(d{k})
                unitWeightVector(l)=1-p;
            else
                unitWeightVector(l)=-(1-p);
%                 unitWeightVector(l)=0;
            end
            l=l+1;
        end
    end
    unitVecs{i}=unitWeightVector;
end
disp('Done with calculating unit weight vectors');
end

function [varargout]=psth_wStd(spikes,binsize,bsmooth,duration,nTrials,theseTrials)

if nargin < 2
    binsize = 50; 
end
if nargin < 3
    bsmooth = 1;
end
% Set duration and number of trials
if ~isempty(nTrials)
    numtrials=nTrials;
elseif isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

nsForStdev=zeros(numtrials,size(n,2));
if ~isempty(theseTrials)
    allTrials=theseTrials;
else
    allTrials=unique(spikes.trials);
end
if length(allTrials)~=numtrials
    if ~isempty(theseTrials)
        allTrials=theseTrials;
    elseif length(spikes.sweeps.trials)==numtrials
        allTrials=spikes.sweeps.trials;
    else
        disp('Needed to fill in trials -- be sure you are using contiguous daq files');
        allTrials=min(unique(spikes.trials)):max(unique(spikes.trials));
    end
end      
for i=1:length(allTrials)
    cspikes=filtspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
if all(isnan(n))
    n = 0;
end
% Compute center of bins
centers = edges + diff(edges(1:2))/2;
% Last point of n contains values falling on edge(end) -- usually zero
if bsmooth
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),0,1);
end

function [varargout]=calcMeanAndStdForUnit(spikes,window,theseTrials)

binsize=1; % in ms
if isempty(theseTrials)
    disp('theseTrials should not be empty in calcMeanAndStdForUnit');
end
numtrials=length(theseTrials);
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges=window(1):binsize:window(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

nsForStdev=zeros(numtrials,size(n,2));
allTrials=theseTrials;
for i=1:length(allTrials)
    cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
s=std(mean(nsForStdev,2),0,1);
if all(isnan(n))
    n = 0;
end
varargout{1} = mean(mean(nsForStdev,2));
varargout{2} = s;
varargout{3} = mean(nsForStdev,2);
end
