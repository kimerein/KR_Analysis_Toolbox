function [prefAv,nonprefAv,trialByTrial_pref]=prefNonprefRealSpikerate(spikes,psth,useAssigns,orderCells)

getn=8;
usel=[0];
ds=10;

led=psth.unitLED{1};
stimcond=psth.unitStimcond{1};

useAssigns=1:length(useAssigns);

prefAv=zeros(size(orderCells,2),size(psth.psths{1},2));
trialByTrial_pref=cell(1,size(orderCells,2));
for i=1:size(orderCells,2)
    [~,sortInd]=sort(orderCells(:,i));
    useUnits=useAssigns(sortInd(end-getn+1:end));
    sumResponse=zeros(length(useUnits),size(psth.psths{1},2));
    sumTrials=zeros(sum(ismember(led,usel) & ismember(stimcond,i)),size(psth.psths{1},2));
    for j=1:length(useUnits)
        currUnit=useUnits(j);
        currpsth=psth.psths{currUnit};
        avpsth=nanmean(currpsth(ismember(led,usel) & ismember(stimcond,i),:),1);
        sumTrials=sumTrials+currpsth(ismember(led,usel) & ismember(stimcond,i),:);
        sumResponse(j,:)=avpsth;
    end
    prefAv(i,:)=nanmean(sumResponse,1);
    trialByTrial_pref{i}=sumTrials./length(useUnits);
end

nonprefAv=zeros(size(orderCells,2),size(psth.psths{1},2));
for i=1:size(orderCells,2)
    [~,sortInd]=sort(-orderCells(:,i));
    useUnits=useAssigns(sortInd(end-getn+1:end));
    sumResponse=zeros(length(useUnits),size(psth.psths{1},2));
    for j=1:length(useUnits)
        currUnit=useUnits(j);
        currpsth=psth.psths{currUnit};
        avpsth=nanmean(currpsth(ismember(led,usel)& ismember(stimcond,i),:),1);
        sumResponse(j,:)=avpsth;
    end
    nonprefAv(i,:)=nanmean(sumResponse,1);
end
    
figure(); 
plot(downSampAv(psth.t,ds),downSampAv(nanmean(prefAv,1),ds),'Color','r');
hold on;
plot(downSampAv(psth.t,ds),downSampAv(nanmean(nonprefAv,1),ds),'Color','b');

 
    
    
    
    