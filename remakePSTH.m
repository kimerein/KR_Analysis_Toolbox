function [FS_exptData,RS_exptData]=remakePSTH(bestSpikes,newassignsinfo,exptData,saveDir)

FSthresh=0.00032;

% origassigns=[];
% uniTrodes=unique(newassignsinfo.trode);
% a=unique(bestSpikes.assigns);
% for i=1:length(uniTrodes)
%     currTrode=uniTrodes(i);
%     if i==1
%         origassigns=a(newassignsinfo.trode==currTrode);
%     else
%         origassigns=[origassigns a(newassignsinfo.trode==currTrode)-origassigns(end)];
%     end
% end

exptDir='W:\Expt Backups\Experiments\';

uniTrodes=unique(newassignsinfo.trode);

% Load expt and get real channels for these trodes
betterAssignsInfo=newassignsinfo;
cload=load([exptDir exptData.name '_expt.mat']);
expt=cload.expt; 
for i=1:length(uniTrodes)
    cUniTrode=uniTrodes(i);
    thisChs=expt.sort.trode(cUniTrode).channels;
    if thisChs(1)==23
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=1;
    elseif thisChs(1)==19
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=2;
    elseif thisChs(1)==17
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=3;
    elseif thisChs(1)==18
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=4;
    else
        disp('problem with trode assigns');
    end
end

concatassigns=[];
newass_assigns=[];
uniTrodes=unique(betterAssignsInfo.trode);
maxSoFar=0;
for i=1:length(uniTrodes)
    currTrode=uniTrodes(i);
    if i==1
        concatassigns=betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode);
        newass_assigns=betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode);
    else
        maxSoFar=max(concatassigns);
        concatassigns=[concatassigns betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode)+maxSoFar];
        newass_assigns=[newass_assigns betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode)];
    end
end
disp('Should be the same numbers');
disp([unique(bestSpikes.assigns)' concatassigns']);

% FSassigns=newassignsinfo.original_assigns(newassignsinfo.waveformWidths<FSthresh);
% RSassigns=newassignsinfo.original_assigns(newassignsinfo.waveformWidths>=FSthresh);
if ~isfield(newassignsinfo,'isFS')
    FSassigns=[];
    RSassigns=[];
    for i=1:length(concatassigns)
        getnewass_assigns=newass_assigns(i);
        getInd=find(newassignsinfo.original_assigns==getnewass_assigns);
        wavewidth=newassignsinfo.waveformWidths(getInd);
        if wavewidth<FSthresh
            FSassigns=[FSassigns concatassigns(i)];
        else
            RSassigns=[RSassigns concatassigns(i)];
        end
    end
else
    FSassigns=concatassigns(newassignsinfo.isFS==1);
    RSassigns=concatassigns(newassignsinfo.isFS==0);
    disp([concatassigns; newassignsinfo.original_assigns]);
end    

[x,y1,y2]=compareMUA(filtspikes(bestSpikes,0,'assigns',FSassigns),exptData.AIs,[],[],exptData);
% [x,y1,y2]=scriptForComparingMUA(filtspikes(bestSpikes,0,'assigns',FSassigns),[],[],[],exptData);
FS_exptData=exptData;
FS_exptData.xpoints=x;
FS_exptData.ypoints1=y1;
FS_exptData.ypoints2=y2;
[x,y1,y2]=compareMUA(filtspikes(bestSpikes,0,'assigns',RSassigns),exptData.AIs,[],[],exptData);
RS_exptData=exptData;
RS_exptData.xpoints=x;
RS_exptData.ypoints1=y1;
RS_exptData.ypoints2=y2;
save([saveDir '\RSFS_byeye.mat'],'FS_exptData','RS_exptData');
end

function [xpoints1,ypoints1,ypoints2]=compareMUA(spikes,useTheseFileInds,useBlackTrials,useRedTrials,exptData)

ledValue=exptData.useLEDcond{2}; % red
noLedValues=exptData.useLEDcond{1}; % black
binsize=1; % ms

stimconds=exptData.useStimcond;
nunits=1;

spikes=filtspikes(spikes,0,'stimcond',stimconds); 

if ~isempty(useTheseFileInds)
    temp=[];
    temp1=[];
    for i=1:length(ledValue)
        spikes=makeTempField(spikes,'led',ledValue(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_withLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
    temp=[];
    temp1=[];
    for i=1:length(noLedValues)
        spikes=makeTempField(spikes,'led',noLedValues(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_noLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
else
    temp=[];
    temp1=[];
    for i=1:length(ledValue)
        spikes=makeTempField(spikes,'led',ledValue(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_withLED=filtspikes(spikes,0,'temp',1);
    temp=[];
    temp1=[];
    for i=1:length(noLedValues)
        spikes=makeTempField(spikes,'led',noLedValues(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_noLED=filtspikes(spikes,0,'temp',1);
end

allSpiketimes_withLED=[];
allSpiketimes_noLED=[];

if ~isempty(useBlackTrials)
    allUnits_noLED=filtspikes(allUnits_noLED,0,'trials',sort(useBlackTrials));
end
if ~isempty(useRedTrials)
    allUnits_withLED=filtspikes(allUnits_withLED,0,'trials',sort(useRedTrials));
end
    
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(allUnits_noLED,binsize);

[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psth_wStdev_valuesOnly(allUnits_withLED,binsize);

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;

% figure();
% plot(xpoints1,ypoints1,'Color','black');
% hold on;
% plot(xpoints2,ypoints2,'Color','red');
% hold off;
end