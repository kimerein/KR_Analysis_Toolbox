function [raw,norm,units,histData,waveformWidth,unitWaveforms]=putSilencingDataTogether(useDir,listing)
%listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');

subtractBaseForUnits=1;
weightAverage=1;
plotThalDepAndInd=1;
removePhotoartifact=1;
nHistBins=5;
colorByDepth=0;
excludeEverywhere=1;
exclude=[5 10 12]-2; 
% exclude=[8 16 8-2 6-2 9-2 11-2 14-2 15-2 17-2]; % Set A
% exclude=[8 16 8-2 3-2 4-2 5-2 7-2 10-2 12-2 13-2 16-2 18-2]; % Set B
doSpont=0; % let zero be zero
timeBefore=0.5; % time before LED onset to take
timeAfter=0.5; % time after LED onset to take
% baseRelativeToLEDonset=[-1 -0.5]; 
% baseRelativeToLEDonset=[-0.2 -0.05]; 
% baseRelativeToLEDonset=[0 1]; 
baseRelativeToLEDonset=[-0.5 -0.3];
% peakRelativeToLEDonset=[-0.3 0];
peakRelativeToLEDonset=[-0.1 0.2];
downSampFactor=1;
alpha1=0.05;
alpha2=0.01;

cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};

% listing=dir(useDir);
% Put PSTH data together
indBefore=timeBefore*1000;
indAfter=timeAfter*1000;
useBase=baseRelativeToLEDonset.*1000;
usePeak=peakRelativeToLEDonset.*1000;
% useBase=baseRelativeToLEDonset.*1000+timeBefore*1000;
% usePeak=peakRelativeToLEDonset.*1000+timeBefore*1000;
ncon=zeros(1,length(listing));
nled=zeros(1,length(listing));
allcon=zeros(length(listing),indBefore+indAfter);
allled=zeros(length(listing),indBefore+indAfter);
useX=-indBefore+1:1:indAfter; % in ms
for i=1:length(listing)
    currFolder=listing(i).name;
    a=load([useDir '\' currFolder '\exptData.mat']);
    exptData=a.exptData;
    ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
    firstInd=find(exptData.xpoints>=ledOnset,1);
    curry1=exptData.ypoints1;
    curry2=exptData.ypoints2;
    if removePhotoartifact==1
        curry1(firstInd-1:firstInd)=ones(1,2).*curry1(firstInd-2);
        curry2(firstInd-1:firstInd)=ones(1,2).*curry2(firstInd-2);
    end
    allcon(i,:)=curry1(firstInd-indBefore+1:firstInd+indAfter);
    allled(i,:)=curry2(firstInd-indBefore+1:firstInd+indAfter);
    ncon(i)=exptData.numtrials1;
    nled(i)=exptData.numtrials2;
end

% Make weighted average
scon=zeros(1,indBefore+indAfter);
sled=zeros(1,indBefore+indAfter);
nusing=1;
for i=1:length(ncon)
    if ismember(i,exclude)
        disp('excluding');
        disp(listing(i).name);
        continue
    end
    if weightAverage==1
        scon=scon+ncon(i)*allcon(i,:);
        sled=sled+nled(i)*allled(i,:);
    else
        scon=scon+allcon(i,:);
        sled=sled+allled(i,:);
        nusing=nusing+1;
    end
end
if weightAverage==1
    scon=scon/sum(ncon);
    sled=sled/sum(nled);
else
    scon=scon/nusing;
    sled=sled/nusing;
end
raw.data1=scon;
raw.data2=sled;
raw.x=useX-timeBefore;
figure(); 
plot(useX-timeBefore,raw.data1,'Color','k');
hold on;
plot(useX-timeBefore,raw.data2,'Color','r');

% Make normalized average
if doSpont==1
    [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
else
    [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
end
norm.x=newX;
norm.y1=newy1;
norm.y1_std=s1;
norm.y2=newy2;
norm.y2_std=s2;

% Put unit data together
led1UnitData=[];
led2UnitData=[];
baseUnitData1=[];
baseUnitData2=[];
calibratedEvCh=[];
waveformWidth=[];
setOfWaveforms=zeros(1000,48);
k=1;
for i=1:length(listing)
    currFolder=listing(i).name;
    if exist([useDir '\' currFolder '\unitsFx.mat'],'file')
        a=load([useDir '\' currFolder '\unitsFx.mat']);
        b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
        if subtractBaseForUnits==1
            c=load([useDir '\' currFolder '\baseUnitsFx.mat']);
        end
    else
        continue
    end
    if excludeEverywhere==1 && ismember(i,exclude)
        continue
    end
    unitsFx=a.unitsFx;
    if subtractBaseForUnits==1
        baseUnits=c.unitsFx;
    end
    newassignsinfo=b.newassignsinfo;
    led1UnitData=[led1UnitData unitsFx.led1Spikes];
    led2UnitData=[led2UnitData unitsFx.led2Spikes];
    if subtractBaseForUnits==1
        baseUnitData1=[baseUnitData1 baseUnits.led1Spikes];
        baseUnitData2=[baseUnitData2 baseUnits.led2Spikes];
    end
    calibratedEvCh=[calibratedEvCh newassignsinfo.calibrated_evCh];
    waveformWidth=[waveformWidth newassignsinfo.waveformWidths];
    setOfWaveforms(k:k+size(newassignsinfo.waveforms,1)-1,:)=newassignsinfo.waveforms(:,1:48);
    k=k+size(newassignsinfo.waveforms,1);
end
unitWaveforms=setOfWaveforms(1:k-1,:);
if subtractBaseForUnits==1
    led1UnitData=led1UnitData-baseUnitData1;
    led2UnitData=led2UnitData-baseUnitData2;
end
units.led1=led1UnitData;
units.led2=led2UnitData;
figure();
if colorByDepth==1
    k=1;
    for i=1:4:16
        scatter(led1UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),led2UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
        hold on;
        k=k+1;
    end
else
    scatter(led1UnitData,led2UnitData);
end
xlabel('LED off -- Hz');
ylabel('LED on -- Hz');

thalDep=led1UnitData-led2UnitData;
thalIndep=led2UnitData;
if plotThalDepAndInd==1
    figure();
    if colorByDepth==1
        k=1;
        for i=1:4:16
            scatter(thalDep(calibratedEvCh>=i & calibratedEvCh<=i+3),thalIndep(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
            hold on;
            k=k+1;
        end
    else
        scatter(thalDep,thalIndep);
    end
    xlabel('Thal. Dependent -- Hz');
    ylabel('Thal. Independent -- Hz');
end

% Make histogram of fractional suppresion
fracSupps=[];
for i=1:length(listing)
    if ismember(i,exclude)
        continue
    end
    currFolder=listing(i).name;
    if exist([useDir '\' currFolder '\fracSupp.mat'],'file')
        a=load([useDir '\' currFolder '\fracSupp.mat']);
    else
        continue
    end
    fracSupp=a.fracSupp;
    fracSupps=[fracSupps fracSupp];
end
[n,xout]=hist(fracSupps,nHistBins);
figure(); 
plot(xout,n);
histData.xout=xout;
histData.n=n;

end

function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putTogetherSilencingPSTH(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2,keepZero,base,peak,scale)

% Use a larger bin size for PSTHs
for i=1:size(dataY1,1)
    [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
end

% Normalize PSTHs
baselineInds=find(newXpoints>base(1) & newXpoints<base(2));
peakInds=find(newXpoints>peak(1) & newXpoints<peak(2));
sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% Normalize by trial WITHOUT LED
for i=1:size(newDataY1,1)
    if keepZero==1
        height=mean(newDataY1(i,baselineInds));
    else
        height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
    end
    if scale==1
        normFactor(i)=1/height;
    else
        normFactor(i)=1;
    end
end
for i=1:size(newDataY1,1)
    if keepZero==1
        newDataY1(i,:)=newDataY1(i,:)*normFactor(i);
        newDataY2(i,:)=newDataY2(i,:)*normFactor(i);
    else
        newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
        newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
    end
end
% newDataY1=newDataY1+sharedBaseline;
% newDataY2=newDataY2+sharedBaseline;

% Get std. dev. measurements
stds1=std(newDataY1,[],1);
stds2=std(newDataY2,[],1);

% All lines figures
figure();
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0.15 0 0.7],[0.7 0.95 0.8],[0.95 0.7 0.8],[0.8 0.65 0.7],[0.35 0.7 0.8],[0.2 0.45 0.2]...
    [0.15 0.4 0.73],[0.3 0.95 0.83],[0.95 0.3 0.83],[0.8 0.35 0.73],[0.35 0.2 0.83],[0.2 0.25 0.23]};
% color legend
for i=1:length(cs)
    line([i i],[0 1],'Color',cs{i});
end
figure();
for i=1:size(newDataY1,1)
    plot(newXpoints,newDataY1(i,:),'Color',cs{i});
    hold on;
end
figure();
for i=1:size(newDataY2,1)
    plot(newXpoints,newDataY2(i,:),'Color',cs{i});
    hold on;
end

figure();
hax=axes();
hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
hold on;
hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);

pvals=zeros(length(newXpoints));
for i=1:length(newXpoints)
    [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
    pvals(i)=p;
end

maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
halfX=(newXpoints(2)-newXpoints(1))/2;
for i=1:length(newXpoints)
    if pvals(i)<alpha1
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
    end
    above=(maxBoth-minBoth)/15;
    if pvals(i)<alpha2
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
    end
end
end



% function [raw,norm,units,histData,waveformWidth,unitWaveforms]=putSilencingDataTogether(useDir,listing)
% %listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');
% 
% subtractBaseForUnits=1;
% weightAverage=0;
% plotThalDepAndInd=0;
% nHistBins=4;
% colorByDepth=0;
% excludeEverywhere=1;
% exclude=[]; %15
% % useWeights=[1 2 3 3 3 3 3];
% useWeights=[];
% % exclude=[8 16 8-2 6-2 9-2 11-2 14-2 15-2 17-2]; % Set A
% % exclude=[8 16 8-2 3-2 4-2 5-2 7-2 10-2 12-2 13-2 16-2 18-2]; % Set B
% doSpont=0; % let zero be zero
% timeBefore=0.5; % time before LED onset to take
% timeAfter=0.5; % time after LED onset to take
% % baseRelativeToLEDonset=[-1 -0.5]; 
% % baseRelativeToLEDonset=[-0.2 -0.05]; 
% % baseRelativeToLEDonset=[0 0.5]; 
% baseRelativeToLEDonset=[-0.5 -0.25];
% % peakRelativeToLEDonset=[-0.3 0];
% peakRelativeToLEDonset=[-0.2 -0.1];
% downSampFactor=1;
% % downSampFactor=22;
% downSampRaw=0;
% alpha1=0.05;
% alpha2=0.01;
% 
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% 
% % listing=dir(useDir);
% % Put PSTH data together
% indBefore=timeBefore*1000;
% indAfter=timeAfter*1000;
% useBase=baseRelativeToLEDonset.*1000;
% usePeak=peakRelativeToLEDonset.*1000;
% % useBase=baseRelativeToLEDonset.*1000+timeBefore*1000;
% % usePeak=peakRelativeToLEDonset.*1000+timeBefore*1000;
% ncon=zeros(1,length(listing));
% nled=zeros(1,length(listing));
% allcon=zeros(length(listing),indBefore+indAfter);
% allled=zeros(length(listing),indBefore+indAfter);
% useX=-indBefore+1:1:indAfter; % in ms
% for i=1:length(listing)
%     currFolder=listing(i).name;
%     if ~exist([useDir '\' currFolder '\exptData.mat'],'file')
%         continue
%     end
%     a=load([useDir '\' currFolder '\exptData.mat']);
% %     if ~exist([useDir '\' currFolder '\FSRS.mat'],'file')
% %         continue
% %     end
% %     a=load([useDir '\' currFolder '\FSRS.mat']);
%     exptData=a.exptData;
% %     exptData=a.RS_exptData;
%     ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
%     firstInd=find(exptData.xpoints>=ledOnset,1);
%     allcon(i,:)=exptData.ypoints1(firstInd-indBefore+1:firstInd+indAfter);
%     allled(i,:)=exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
%     ncon(i)=exptData.numtrials1;
%     nled(i)=exptData.numtrials2;
% end
% 
% % Make weighted average
% scon=zeros(1,indBefore+indAfter);
% sled=zeros(1,indBefore+indAfter);
% nusing=1;
% for i=1:length(ncon)
%     if ismember(i,exclude)
%         disp('excluding');
%         disp(listing(i).name);
%         continue
%     end
%     if weightAverage==1 && isempty(useWeights)
%         scon=scon+ncon(i)*allcon(i,:);
%         sled=sled+nled(i)*allled(i,:);
%     elseif weightAverage==1
%         scon=scon+useWeights(i)*allcon(i,:);
%         sled=sled+useWeights(i)*allled(i,:);
%     else
%         scon=scon+allcon(i,:);
%         sled=sled+allled(i,:);
%         nusing=nusing+1;
%     end
% end
% if weightAverage==1
%     scon=scon/sum(ncon);
%     sled=sled/sum(nled);
% else
%     scon=scon/nusing;
%     sled=sled/nusing;
% end
% if downSampRaw==1
%     raw.data1=downSampAv(scon,downSampFactor);
%     raw.data2=downSampAv(sled,downSampFactor);
%     raw.x=downSampAv(useX,downSampFactor)-timeBefore;
% else
%     raw.data1=scon;
%     raw.data2=sled;
%     raw.x=useX-timeBefore;
% end
% figure(); 
% plot(raw.x,raw.data1,'Color','k');
% hold on;
% plot(raw.x,raw.data2,'Color','r');
% 
% % Make normalized average
% if doSpont==1
%     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% else
%     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% end
% norm.x=newX;
% norm.y1=newy1;
% norm.y1_std=s1;
% norm.y2=newy2;
% norm.y2_std=s2;
% 
% % Put unit data together
% led1UnitData=[];
% led2UnitData=[];
% baseUnitData1=[];
% baseUnitData2=[];
% calibratedEvCh=[];
% waveformWidth=[];
% setOfWaveforms=zeros(1000,48);
% k=1;
% for i=1:length(listing)
%     currFolder=listing(i).name;
%     if exist([useDir '\' currFolder '\unitsFx.mat'],'file')
%         a=load([useDir '\' currFolder '\unitsFx.mat']);
%         b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
%         if subtractBaseForUnits==1
%             c=load([useDir '\' currFolder '\baseUnitsFx.mat']);
%         end
%     else
%         continue
%     end
%     if excludeEverywhere==1 && ismember(i,exclude)
%         continue
%     end
%     unitsFx=a.unitsFx;
%     if subtractBaseForUnits==1
%         baseUnits=c.unitsFx;
%     end
%     newassignsinfo=b.newassignsinfo;
%     led1UnitData=[led1UnitData unitsFx.led1Spikes];
%     led2UnitData=[led2UnitData unitsFx.led2Spikes];
%     if subtractBaseForUnits==1
%         baseUnitData1=[baseUnitData1 baseUnits.led1Spikes];
%         baseUnitData2=[baseUnitData2 baseUnits.led2Spikes];
%     end
%     calibratedEvCh=[calibratedEvCh newassignsinfo.calibrated_evCh];
%     waveformWidth=[waveformWidth newassignsinfo.waveformWidths];
%     setOfWaveforms(k:k+size(newassignsinfo.waveforms,1)-1,:)=newassignsinfo.waveforms(:,1:48);
%     k=k+size(newassignsinfo.waveforms,1);
% end
% unitWaveforms=setOfWaveforms(1:k-1,:);
% if subtractBaseForUnits==1
%     led1UnitData=led1UnitData-baseUnitData1;
%     led2UnitData=led2UnitData-baseUnitData2;
% end
% units.led1=led1UnitData;
% units.led2=led2UnitData;
% figure();
% if colorByDepth==1
%     k=1;
%     for i=1:4:16
%         scatter(led1UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),led2UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
%         hold on;
%         k=k+1;
%     end
% else
%     scatter(led1UnitData,led2UnitData);
% end
% xlabel('LED off -- Hz');
% ylabel('LED on -- Hz');
% 
% thalDep=led1UnitData-led2UnitData;
% thalIndep=led2UnitData;
% if plotThalDepAndInd==1
%     figure();
%     if colorByDepth==1
%         k=1;
%         for i=1:4:16
%             scatter(thalDep(calibratedEvCh>=i & calibratedEvCh<=i+3),thalIndep(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
%             hold on;
%             k=k+1;
%         end
%     else
%         scatter(thalDep,thalIndep);
%     end
%     xlabel('Thal. Dependent -- Hz');
%     ylabel('Thal. Independent -- Hz');
% end
% 
% % Make histogram of fractional suppresion
% fracSupps=[];
% for i=1:length(listing)
%     if ismember(i,exclude)
%         continue
%     end
%     currFolder=listing(i).name;
%     if exist([useDir '\' currFolder '\fracSupp.mat'],'file')
%         a=load([useDir '\' currFolder '\fracSupp.mat']);
%     else
%         continue
%     end
%     fracSupp=a.fracSupp;
%     fracSupps=[fracSupps fracSupp];
% end
% [n,xout]=hist(fracSupps,nHistBins);
% figure(); 
% plot(xout,n);
% histData.xout=xout;
% histData.n=n;
% 
% end
% 
% function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putTogetherSilencingPSTH(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2,keepZero,base,peak,scale)
% 
% % Use a larger bin size for PSTHs
% for i=1:size(dataY1,1)
%     [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
% end
% 
% % Normalize PSTHs
% baselineInds=find(newXpoints>base(1) & newXpoints<base(2));
% peakInds=find(newXpoints>peak(1) & newXpoints<peak(2));
% sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% % Normalize by trial WITHOUT LED
% for i=1:size(newDataY1,1)
%     if keepZero==1
%         height=mean(newDataY1(i,baselineInds));
%     else
%         height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
%     end
%     if scale==1
%         normFactor(i)=1/height;
%     else
%         normFactor(i)=1;
%     end
% end
% for i=1:size(newDataY1,1)
%     if keepZero==1
%         newDataY1(i,:)=newDataY1(i,:)*normFactor(i);
%         newDataY2(i,:)=newDataY2(i,:)*normFactor(i);
%     else
%         newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
%         newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
%     end
% end
% % newDataY1=newDataY1+sharedBaseline;
% % newDataY2=newDataY2+sharedBaseline;
% 
% % Get std. dev. measurements
% stds1=std(newDataY1,[],1);
% stds2=std(newDataY2,[],1);
% 
% % All lines figures
% figure();
% % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% %     [0.15 0 0.7],[0.7 0.95 0.8],[0.95 0.7 0.8],[0.8 0.65 0.7],[0.35 0.7 0.8],[0.2 0.45 0.2]};
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% 
% % color legend
% for i=1:length(cs)
%     line([i i],[0 1],'Color',cs{i});
% end
% figure();
% for i=1:size(newDataY1,1)
%     plot(newXpoints,newDataY1(i,:),'Color',cs{i});
%     hold on;
% end
% figure();
% for i=1:size(newDataY2,1)
%     plot(newXpoints,newDataY2(i,:),'Color',cs{i});
%     hold on;
% end
% 
% figure();
% hax=axes();
% hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
% addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
% hold on;
% hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
% addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);
% 
% pvals=zeros(length(newXpoints));
% for i=1:length(newXpoints)
%     [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
%     pvals(i)=p;
% end
% 
% maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
% minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
% halfX=(newXpoints(2)-newXpoints(1))/2;
% for i=1:length(newXpoints)
%     if pvals(i)<alpha1
%         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
%     end
%     above=(maxBoth-minBoth)/15;
%     if pvals(i)<alpha2
%         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
%     end
% end
% end
% 
% 
% 
% % function [raw,norm,units,histData,waveformWidth,unitWaveforms,RSvsFS,nonnorm,fracs]=putSilencingDataTogether(useDir,listing)
% % %listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');
% % 
% % raw=[]; norm=[]; units=[]; histData=[]; waveformWidth=[]; unitWaveforms=[];
% % subtractBaseForUnits=1;
% % weightAverage=0;
% % plotThalDepAndInd=0;
% % nHistBins=5;
% % colorByDepth=0;
% % excludeEverywhere=0;
% % % exclude=[9 11 23]; %15
% % exclude=[]; %15
% % % exclude=[8 16 8-2 6-2 9-2 11-2 14-2 15-2 17-2]; % Set A
% % % exclude=[8 16 8-2 3-2 4-2 5-2 7-2 10-2 12-2 13-2 16-2 18-2]; % Set B
% % doSpont=0; % let zero be zero
% % timeBefore=0.5; % time before LED onset to take
% % timeAfter=1; % time after LED onset to take
% % % baseRelativeToLEDonset=[-1 -0.5]; 
% % % baseRelativeToLEDonset=[-0.2 -0.05]; 
% % % baseRelativeToLEDonset=[0 1]; 
% % % baseRelativeToLEDonset=[-1 0];
% % baseRelativeToLEDonset=[-0.5 -0.4];
% % % peakRelativeToLEDonset=[-0.3 0];
% % peakRelativeToLEDonset=[0 0.2];
% % % downSampFactor=12;
% % downSampFactor=20;
% % alpha1=0.05;
% % alpha2=0.01;
% % 
% % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % 
% % % listing=dir(useDir);
% % % Put PSTH data together
% % indBefore=timeBefore*1000;
% % indAfter=timeAfter*1000;
% % useBase=baseRelativeToLEDonset.*1000;
% % usePeak=peakRelativeToLEDonset.*1000;
% % % useBase=baseRelativeToLEDonset.*1000+timeBefore*1000;
% % % usePeak=peakRelativeToLEDonset.*1000+timeBefore*1000;
% % ncon=zeros(1,length(listing));
% % nled=zeros(1,length(listing));
% % allcon=zeros(length(listing),indBefore+indAfter);
% % allled=zeros(length(listing),indBefore+indAfter);
% % FSshutoff=zeros(1,indBefore+indAfter);
% % RSshutoff=zeros(1,indBefore+indAfter);
% % useX=-indBefore+1:1:indAfter; % in ms
% % countingFSRS=0;
% % isFS_together=[];
% % relativeToFS_taus=[];
% % fracs=[];
% % for i=1:length(listing)
% %     currFolder=listing(i).name;
% %     if ~exist([useDir '\' currFolder '\fracSupp.mat'],'file')
% %         continue
% %     end
% %     a=load([useDir '\' currFolder '\fracSupp.mat']);
% %     fracs=[fracs a.fracSupp];
% %     if ~exist([useDir '\' currFolder '\exptData.mat'],'file')
% %         continue
% %     end
% %     a=load([useDir '\' currFolder '\exptData.mat']);
% %     exptData=a.exptData;
% %     ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
% %     firstInd=find(exptData.xpoints>=ledOnset,1);
% %     allcon(i,:)=exptData.ypoints1(firstInd-indBefore+1:firstInd+indAfter);
% %     allled(i,:)=exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
% %     ncon(i)=exptData.numtrials1;
% %     nled(i)=exptData.numtrials2;
% %     if exist([useDir '\' currFolder '\RSFS_byeye.mat'],'file')
% %         a=load([useDir '\' currFolder '\RSFS_byeye.mat']);
% %         b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
% %         nn=b.newassignsinfo;
% %         FS_exptData=a.FS_exptData;
% %         RS_exptData=a.RS_exptData;
% %         firstInd=find(FS_exptData.xpoints>=ledOnset,1);
% % %         FSshutoff=FSshutoff+FS_exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter)*sum(nn.isFS==1);
% % %         RSshutoff=RSshutoff+RS_exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter)*sum(nn.isFS==0);
% %         FSshutoff=FSshutoff+FS_exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
% %         RSshutoff=RSshutoff+RS_exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
% %         countingFSRS=countingFSRS+1;
% %     end
% %     if exist([useDir '\' currFolder '\bestUnitsInfo.mat'],'file') && ...
% %        exist([useDir '\' currFolder '\units_shutoffBETTER.mat'],'file')
% %         a=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
% %         bestUnits_newass=a.newassignsinfo;
% %         a=load([useDir '\' currFolder '\units_shutoffBETTER.mat']);
% %         forFS_allUnitsTimeCourse=a.allUnitsTimeCourse;
% %         forFS_allUnitsTimeCourse.taus=forFS_allUnitsTimeCourse.taus';
% %         takeThese=~isnan(forFS_allUnitsTimeCourse.taus) & ~isnan(bestUnits_newass.isFS);
% %         isFS_together=[isFS_together bestUnits_newass.isFS(takeThese)];
% %         relativeToFS_taus=[relativeToFS_taus forFS_allUnitsTimeCourse.taus(takeThese)];
% %     end
% % end
% % FSshutoff=FSshutoff./countingFSRS;
% % RSshutoff=RSshutoff./countingFSRS;
% % RSvsFS.FSshutoff=FSshutoff;
% % RSvsFS.RSshutoff=RSshutoff;
% % RSvsFS.x=useX-timeBefore;
% % figure(); 
% % plot(useX-timeBefore,RSshutoff);
% % title('RS Average PSTH');
% % figure(); 
% % plot(useX-timeBefore,FSshutoff);
% % title('FS Average PSTH');
% % 
% % figure(); 
% % isFS_together=isFS_together(relativeToFS_taus<0.2);
% % relativeToFS_taus=relativeToFS_taus(relativeToFS_taus<0.2);
% % scatter(isFS_together(relativeToFS_taus<0.2),relativeToFS_taus(relativeToFS_taus<0.2));
% % disp('number of FS units');
% % disp(sum(isFS_together==1));
% % disp('number of RS units');
% % disp(sum(isFS_together==0));
% % title('isFS vs. unit tau');
% % [n,xout]=hist(relativeToFS_taus(isFS_together==1),10);
% % figure(); plot(xout,n,'Color','b');
% % hold on; 
% % [n,xout]=hist(relativeToFS_taus(isFS_together==0),100);
% % plot(xout,n,'Color','r');
% % % return
% % % Make weighted average
% % scon=zeros(1,indBefore+indAfter);
% % sled=zeros(1,indBefore+indAfter);
% % nusing=1;
% % for i=1:length(ncon)
% %     if ismember(i,exclude)
% %         disp('excluding');
% %         disp(listing(i).name);
% %         continue
% %     end
% %     if weightAverage==1
% %         scon=scon+ncon(i)*allcon(i,:);
% %         sled=sled+nled(i)*allled(i,:);
% %     else
% %         scon=scon+allcon(i,:);
% %         sled=sled+allled(i,:);
% %         nusing=nusing+1;
% %     end
% % end
% % if weightAverage==1
% %     scon=scon/sum(ncon);
% %     sled=sled/sum(nled);
% % else
% %     scon=scon/nusing;
% %     sled=sled/nusing;
% % end
% % raw.data1=scon;
% % raw.data2=sled;
% % raw.x=useX-timeBefore;
% % figure(); 
% % plot(useX-timeBefore,raw.data1,'Color','k');
% % hold on;
% % plot(useX-timeBefore,raw.data2,'Color','r');
% % 
% % % Make normalized average
% % nonnorm.y1=allcon(~ismember(1:size(allcon,1),exclude),:);
% % nonnorm.y2=allled(~ismember(1:size(allcon,1),exclude),:);
% % if doSpont==1
% %     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% % else
% %     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% % end
% % norm.x=newX;
% % norm.y1=newy1;
% % norm.y1_std=s1;
% % norm.y2=newy2;
% % norm.y2_std=s2;
% % 
% % % Put unit data together
% % led1UnitData=[];
% % led2UnitData=[];
% % baseUnitData1=[];
% % baseUnitData2=[];
% % calibratedEvCh=[];
% % waveformWidth=[];
% % setOfWaveforms=zeros(1000,48);
% % k=1;
% % for i=1:length(listing)
% %     currFolder=listing(i).name;
% %     if exist([useDir '\' currFolder '\unitsFx.mat'],'file') && exist([useDir '\' currFolder '\bestUnitsInfo.mat'],'file')
% %         a=load([useDir '\' currFolder '\unitsFx.mat']);
% %         b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
% %         if subtractBaseForUnits==1
% %             c=load([useDir '\' currFolder '\baseUnitsFx.mat']);
% %         end
% %     else
% %         continue
% %     end
% %     if excludeEverywhere==1 && ismember(i,exclude)
% %         continue
% %     end
% %     unitsFx=a.unitsFx;
% %     if subtractBaseForUnits==1
% %         baseUnits=c.unitsFx;
% %     end
% %     newassignsinfo=b.newassignsinfo;
% %     led1UnitData=[led1UnitData unitsFx.led1Spikes];
% %     led2UnitData=[led2UnitData unitsFx.led2Spikes];
% %     if subtractBaseForUnits==1
% %         baseUnitData1=[baseUnitData1 baseUnits.led1Spikes];
% %         baseUnitData2=[baseUnitData2 baseUnits.led2Spikes];
% %     end
% %     calibratedEvCh=[calibratedEvCh newassignsinfo.calibrated_evCh];
% %     waveformWidth=[waveformWidth newassignsinfo.waveformWidths];
% %     setOfWaveforms(k:k+size(newassignsinfo.waveforms,1)-1,:)=newassignsinfo.waveforms(:,1:48);
% %     k=k+size(newassignsinfo.waveforms,1);
% % end
% % unitWaveforms=setOfWaveforms(1:k-1,:);
% % if subtractBaseForUnits==1
% %     led1UnitData=led1UnitData-baseUnitData1;
% %     led2UnitData=led2UnitData-baseUnitData2;
% % end
% % units.led1=led1UnitData;
% % units.led2=led2UnitData;
% % figure();
% % if colorByDepth==1
% %     k=1;
% %     for i=1:4:16
% %         scatter(led1UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),led2UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
% %         hold on;
% %         k=k+1;
% %     end
% % else
% %     scatter(led1UnitData,led2UnitData);
% % end
% % xlabel('LED off -- Hz');
% % ylabel('LED on -- Hz');
% % 
% % thalDep=led1UnitData-led2UnitData;
% % thalIndep=led2UnitData;
% % if plotThalDepAndInd==1
% %     figure();
% %     if colorByDepth==1
% %         k=1;
% %         for i=1:4:16
% %             scatter(thalDep(calibratedEvCh>=i & calibratedEvCh<=i+3),thalIndep(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
% %             hold on;
% %             k=k+1;
% %         end
% %     else
% %         scatter(thalDep,thalIndep);
% %     end
% %     xlabel('Thal. Dependent -- Hz');
% %     ylabel('Thal. Independent -- Hz');
% % end
% % 
% % % Make histogram of fractional suppresion
% % fracSupps=[];
% % for i=1:length(listing)
% %     if ismember(i,exclude)
% %         continue
% %     end
% %     currFolder=listing(i).name;
% %     if exist([useDir '\' currFolder '\fracSupp.mat'],'file')
% %         a=load([useDir '\' currFolder '\fracSupp.mat']);
% %     else
% %         continue
% %     end
% %     fracSupp=a.fracSupp;
% %     fracSupps=[fracSupps fracSupp];
% % end
% % [n,xout]=hist(fracSupps,nHistBins);
% % figure(); 
% % plot(xout,n);
% % histData.xout=xout;
% % histData.n=n;
% % 
% % end
% % 
% % function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putTogetherSilencingPSTH(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2,keepZero,base,peak,scale)
% % 
% % % Use a larger bin size for PSTHs
% % for i=1:size(dataY1,1)
% %     [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
% % end
% % 
% % % Normalize PSTHs
% % baselineInds=find(newXpoints>base(1) & newXpoints<base(2));
% % peakInds=find(newXpoints>peak(1) & newXpoints<peak(2));
% % sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% % % Normalize by trial WITHOUT LED
% % for i=1:size(newDataY1,1)
% %     if keepZero==1
% %         height=mean(newDataY1(i,baselineInds));
% %     else
% %         height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
% %     end
% %     if scale==1
% %         normFactor(i)=1/height;
% %     else
% %         normFactor(i)=1;
% %     end
% % end
% % for i=1:size(newDataY1,1)
% %     if keepZero==1
% %         newDataY1(i,:)=newDataY1(i,:)*normFactor(i);
% %         newDataY2(i,:)=newDataY2(i,:)*normFactor(i);
% %     else
% %         newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
% %         newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
% %     end
% % end
% % % newDataY1=newDataY1+sharedBaseline;
% % % newDataY2=newDataY2+sharedBaseline;
% % 
% % % Get std. dev. measurements
% % stds1=std(newDataY1,[],1);
% % stds2=std(newDataY2,[],1);
% % 
% % % All lines figures
% % figure();
% % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% %     [0.15 0 0.7],[0.7 0.95 0.8],[0.95 0.7 0.8],[0.8 0.65 0.7],[0.35 0.7 0.8],[0.2 0.45 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% % % color legend
% % for i=1:length(cs)
% %     line([i i],[0 1],'Color',cs{i});
% % end
% % figure();
% % for i=1:size(newDataY1,1)
% %     plot(newXpoints,newDataY1(i,:),'Color',cs{i});
% %     hold on;
% % end
% % figure();
% % for i=1:size(newDataY2,1)
% %     plot(newXpoints,newDataY2(i,:),'Color',cs{i});
% %     hold on;
% % end
% % 
% % figure();
% % hax=axes();
% % hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
% % addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
% % hold on;
% % hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
% % addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);
% % 
% % pvals=zeros(length(newXpoints));
% % for i=1:length(newXpoints)
% %     [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
% %     pvals(i)=p;
% % end
% % 
% % maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
% % minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
% % halfX=(newXpoints(2)-newXpoints(1))/2;
% % for i=1:length(newXpoints)
% %     if pvals(i)<alpha1
% %         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
% %     end
% %     above=(maxBoth-minBoth)/15;
% %     if pvals(i)<alpha2
% %         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
% %     end
% % end
% % end
% % % function [raw,norm,units,histData,waveformWidth,unitWaveforms]=putSilencingDataTogether(useDir,listing)
% % % %listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');
% % % 
% % % subtractBaseForUnits=0;
% % % weightAverage=0;
% % % plotThalDepAndInd=0;
% % % nHistBins=3;
% % % colorByDepth=0;
% % % excludeEverywhere=1;
% % % exclude=[]; %15
% % % % useWeights=[1 2 3 3 3 3 3];
% % % useWeights=[];
% % % % exclude=[8 16 8-2 6-2 9-2 11-2 14-2 15-2 17-2]; % Set A
% % % % exclude=[8 16 8-2 3-2 4-2 5-2 7-2 10-2 12-2 13-2 16-2 18-2]; % Set B
% % % doSpont=1; % let zero be zero
% % % timeBefore=0.5; % time before LED onset to take
% % % timeAfter=0.5; % time after LED onset to take
% % % % baseRelativeToLEDonset=[-1 -0.5]; 
% % % % baseRelativeToLEDonset=[-0.2 -0.05]; 
% % % baseRelativeToLEDonset=[0 0.5]; 
% % % % baseRelativeToLEDonset=[-0.5 -0.25];
% % % % peakRelativeToLEDonset=[-0.3 0];
% % % peakRelativeToLEDonset=[-0.15 0];
% % % % downSampFactor=30;
% % % downSampFactor=22;
% % % downSampRaw=0;
% % % alpha1=0.05;
% % % alpha2=0.01;
% % % 
% % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % % 
% % % % listing=dir(useDir);
% % % % Put PSTH data together
% % % indBefore=timeBefore*1000;
% % % indAfter=timeAfter*1000;
% % % useBase=baseRelativeToLEDonset.*1000;
% % % usePeak=peakRelativeToLEDonset.*1000;
% % % % useBase=baseRelativeToLEDonset.*1000+timeBefore*1000;
% % % % usePeak=peakRelativeToLEDonset.*1000+timeBefore*1000;
% % % ncon=zeros(1,length(listing));
% % % nled=zeros(1,length(listing));
% % % allcon=zeros(length(listing),indBefore+indAfter);
% % % allled=zeros(length(listing),indBefore+indAfter);
% % % useX=-indBefore+1:1:indAfter; % in ms
% % % for i=1:length(listing)
% % %     currFolder=listing(i).name;
% % %     if ~exist([useDir '\' currFolder '\exptData.mat'],'file')
% % %         continue
% % %     end
% % %     a=load([useDir '\' currFolder '\exptData.mat']);
% % % %     if ~exist([useDir '\' currFolder '\FSRS.mat'],'file')
% % % %         continue
% % % %     end
% % % %     a=load([useDir '\' currFolder '\FSRS.mat']);
% % %     exptData=a.exptData;
% % % %     exptData=a.RS_exptData;
% % %     ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
% % %     firstInd=find(exptData.xpoints>=ledOnset,1);
% % %     allcon(i,:)=exptData.ypoints1(firstInd-indBefore+1:firstInd+indAfter);
% % %     allled(i,:)=exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
% % %     ncon(i)=exptData.numtrials1;
% % %     nled(i)=exptData.numtrials2;
% % % end
% % % 
% % % % Make weighted average
% % % scon=zeros(1,indBefore+indAfter);
% % % sled=zeros(1,indBefore+indAfter);
% % % nusing=1;
% % % for i=1:length(ncon)
% % %     if ismember(i,exclude)
% % %         disp('excluding');
% % %         disp(listing(i).name);
% % %         continue
% % %     end
% % %     if weightAverage==1 && isempty(useWeights)
% % %         scon=scon+ncon(i)*allcon(i,:);
% % %         sled=sled+nled(i)*allled(i,:);
% % %     elseif weightAverage==1
% % %         scon=scon+useWeights(i)*allcon(i,:);
% % %         sled=sled+useWeights(i)*allled(i,:);
% % %     else
% % %         scon=scon+allcon(i,:);
% % %         sled=sled+allled(i,:);
% % %         nusing=nusing+1;
% % %     end
% % % end
% % % if weightAverage==1
% % %     scon=scon/sum(ncon);
% % %     sled=sled/sum(nled);
% % % else
% % %     scon=scon/nusing;
% % %     sled=sled/nusing;
% % % end
% % % if downSampRaw==1
% % %     raw.data1=downSampAv(scon,downSampFactor);
% % %     raw.data2=downSampAv(sled,downSampFactor);
% % %     raw.x=downSampAv(useX,downSampFactor)-timeBefore;
% % % else
% % %     raw.data1=scon;
% % %     raw.data2=sled;
% % %     raw.x=useX-timeBefore;
% % % end
% % % figure(); 
% % % plot(raw.x,raw.data1,'Color','k');
% % % hold on;
% % % plot(raw.x,raw.data2,'Color','r');
% % % 
% % % % Make normalized average
% % % if doSpont==1
% % %     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% % % else
% % %     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% % % end
% % % norm.x=newX;
% % % norm.y1=newy1;
% % % norm.y1_std=s1;
% % % norm.y2=newy2;
% % % norm.y2_std=s2;
% % % 
% % % % Put unit data together
% % % led1UnitData=[];
% % % led2UnitData=[];
% % % baseUnitData1=[];
% % % baseUnitData2=[];
% % % calibratedEvCh=[];
% % % waveformWidth=[];
% % % setOfWaveforms=zeros(1000,48);
% % % k=1;
% % % for i=1:length(listing)
% % %     currFolder=listing(i).name;
% % %     if exist([useDir '\' currFolder '\unitsFx.mat'],'file')
% % %         a=load([useDir '\' currFolder '\unitsFx.mat']);
% % %         b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
% % %         if subtractBaseForUnits==1
% % %             c=load([useDir '\' currFolder '\baseUnitsFx.mat']);
% % %         end
% % %     else
% % %         continue
% % %     end
% % %     if excludeEverywhere==1 && ismember(i,exclude)
% % %         continue
% % %     end
% % %     unitsFx=a.unitsFx;
% % %     if subtractBaseForUnits==1
% % %         baseUnits=c.unitsFx;
% % %     end
% % %     newassignsinfo=b.newassignsinfo;
% % %     led1UnitData=[led1UnitData unitsFx.led1Spikes];
% % %     led2UnitData=[led2UnitData unitsFx.led2Spikes];
% % %     if subtractBaseForUnits==1
% % %         baseUnitData1=[baseUnitData1 baseUnits.led1Spikes];
% % %         baseUnitData2=[baseUnitData2 baseUnits.led2Spikes];
% % %     end
% % %     calibratedEvCh=[calibratedEvCh newassignsinfo.calibrated_evCh];
% % %     waveformWidth=[waveformWidth newassignsinfo.waveformWidths];
% % %     setOfWaveforms(k:k+size(newassignsinfo.waveforms,1)-1,:)=newassignsinfo.waveforms(:,1:48);
% % %     k=k+size(newassignsinfo.waveforms,1);
% % % end
% % % unitWaveforms=setOfWaveforms(1:k-1,:);
% % % if subtractBaseForUnits==1
% % %     led1UnitData=led1UnitData-baseUnitData1;
% % %     led2UnitData=led2UnitData-baseUnitData2;
% % % end
% % % units.led1=led1UnitData;
% % % units.led2=led2UnitData;
% % % figure();
% % % if colorByDepth==1
% % %     k=1;
% % %     for i=1:4:16
% % %         scatter(led1UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),led2UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
% % %         hold on;
% % %         k=k+1;
% % %     end
% % % else
% % %     scatter(led1UnitData,led2UnitData);
% % % end
% % % xlabel('LED off -- Hz');
% % % ylabel('LED on -- Hz');
% % % 
% % % thalDep=led1UnitData-led2UnitData;
% % % thalIndep=led2UnitData;
% % % if plotThalDepAndInd==1
% % %     figure();
% % %     if colorByDepth==1
% % %         k=1;
% % %         for i=1:4:16
% % %             scatter(thalDep(calibratedEvCh>=i & calibratedEvCh<=i+3),thalIndep(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
% % %             hold on;
% % %             k=k+1;
% % %         end
% % %     else
% % %         scatter(thalDep,thalIndep);
% % %     end
% % %     xlabel('Thal. Dependent -- Hz');
% % %     ylabel('Thal. Independent -- Hz');
% % % end
% % % 
% % % % Make histogram of fractional suppresion
% % % fracSupps=[];
% % % for i=1:length(listing)
% % %     if ismember(i,exclude)
% % %         continue
% % %     end
% % %     currFolder=listing(i).name;
% % %     if exist([useDir '\' currFolder '\fracSupp.mat'],'file')
% % %         a=load([useDir '\' currFolder '\fracSupp.mat']);
% % %     else
% % %         continue
% % %     end
% % %     fracSupp=a.fracSupp;
% % %     fracSupps=[fracSupps fracSupp];
% % % end
% % % [n,xout]=hist(fracSupps,nHistBins);
% % % figure(); 
% % % plot(xout,n);
% % % histData.xout=xout;
% % % histData.n=n;
% % % 
% % % end
% % % 
% % % function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putTogetherSilencingPSTH(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2,keepZero,base,peak,scale)
% % % 
% % % % Use a larger bin size for PSTHs
% % % for i=1:size(dataY1,1)
% % %     [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
% % % end
% % % 
% % % % Normalize PSTHs
% % % baselineInds=find(newXpoints>base(1) & newXpoints<base(2));
% % % peakInds=find(newXpoints>peak(1) & newXpoints<peak(2));
% % % sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% % % % Normalize by trial WITHOUT LED
% % % for i=1:size(newDataY1,1)
% % %     if keepZero==1
% % %         height=mean(newDataY1(i,baselineInds));
% % %     else
% % %         height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
% % %     end
% % %     if scale==1
% % %         normFactor(i)=1/height;
% % %     else
% % %         normFactor(i)=1;
% % %     end
% % % end
% % % for i=1:size(newDataY1,1)
% % %     if keepZero==1
% % %         newDataY1(i,:)=newDataY1(i,:)*normFactor(i);
% % %         newDataY2(i,:)=newDataY2(i,:)*normFactor(i);
% % %     else
% % %         newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
% % %         newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
% % %     end
% % % end
% % % % newDataY1=newDataY1+sharedBaseline;
% % % % newDataY2=newDataY2+sharedBaseline;
% % % 
% % % % Get std. dev. measurements
% % % stds1=std(newDataY1,[],1);
% % % stds2=std(newDataY2,[],1);
% % % 
% % % % All lines figures
% % % figure();
% % % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % % %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % % %     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% % % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % % %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% % % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % % %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% % % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % % %     [0.15 0 0.7],[0.7 0.95 0.8],[0.95 0.7 0.8],[0.8 0.65 0.7],[0.35 0.7 0.8],[0.2 0.45 0.2]};
% % % cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
% % %     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
% % %     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% % % 
% % % % color legend
% % % for i=1:length(cs)
% % %     line([i i],[0 1],'Color',cs{i});
% % % end
% % % figure();
% % % for i=1:size(newDataY1,1)
% % %     plot(newXpoints,newDataY1(i,:),'Color',cs{i});
% % %     hold on;
% % % end
% % % figure();
% % % for i=1:size(newDataY2,1)
% % %     plot(newXpoints,newDataY2(i,:),'Color',cs{i});
% % %     hold on;
% % % end
% % % 
% % % figure();
% % % hax=axes();
% % % hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
% % % addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
% % % hold on;
% % % hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
% % % addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);
% % % 
% % % pvals=zeros(length(newXpoints));
% % % for i=1:length(newXpoints)
% % %     [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
% % %     pvals(i)=p;
% % % end
% % % 
% % % maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
% % % minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
% % % halfX=(newXpoints(2)-newXpoints(1))/2;
% % % for i=1:length(newXpoints)
% % %     if pvals(i)<alpha1
% % %         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
% % %     end
% % %     above=(maxBoth-minBoth)/15;
% % %     if pvals(i)<alpha2
% % %         line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
% % %     end
% % % end
% % % end
