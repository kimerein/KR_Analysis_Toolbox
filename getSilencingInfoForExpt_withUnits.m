function getSilencingInfoForExpt_withUnits(expt,orderedByDepthSpikes)
% Start with non-filtered spikes, so can read event_channel out of struct
% without having to recalculate

justExcludeAssigns=0;
doMU=0;

% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN stim. duration 0 ms\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\AWAKE\Evoked LP\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Awake\Drifting Evoked\';
% saveDir='W:\Analysis Computer\Across Mice Finals\V2L\Spont Awake\';
% saveDir='W:\Analysis Computer\Across Mice Finals\V2L\Evoked Awake\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont LP\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Evoked Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Spont Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Static Silencing Across Mice\';
saveDir='W:\Analysis Computer\New Fig 1 150115\Anesth 10 ms Flash Evoked\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Different LED Intensities\';
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
Fs=32000;
spontExpt=0; % 0 if evoked silencing experiment
dLGNdiscard=0; % don't use units with half-width <0.22 ms in frac. supp. and PSTH

Trode_numbers=[2 3 4]; % Must be in order of depth
Trode_orderByDepth=[1 2 3];
 
criterionAssigns{1}=[56]; % This matches Trode_number
criterionAssigns{2}=[1 3 16 22 23 25 26 28 31 32 88 89 90 91 92 97 98 101 102 105 109];
criterionAssigns{3}=[5 7 8 9 10 11 19 21 23 24 25 27 28 29 36 37 38 40 47 100 104 107 118 119 121];
% criterionAssigns{4}=[5];

bestAssigns{1}=[56]; % This matches Trode_number
bestAssigns{2}=[1 3 16 22 23 25 26 28 31 32 88 89 90 91 92 97 98 101 102 105 109];
bestAssigns{3}=[5 7 8 9 10 11 19 21 23 24 25 27 28 29 36 37 38 40 47 100 104 107 118 119 121];
% bestAssigns{4}=[5];
% bestAssigns{1}=[]; % This matches Trode_number
% bestAssigns{2}=[];
% bestAssigns{3}=[];
% bestAssigns{4}=[];

checkLED=0;
aliLED=0;
params.useFileInd=[21:29];
params.useStimcond=[1];
params.useTrigger=1:12;
params.useLEDcond={[5]; [5.05]};
params.trialDuration=5; % in s
params.stimulusOn=[1 3]; % in s relative to trial onset
params.stimulusType='10X Pulse Train';
params.stimulusDuration=-500; % in ms, stimulus duration before LED onset
params.LEDduration=2500; % in ms, duration of LED pulse -- maybe short to exclude photo-artifact
params.stimulusContrast=1;
params.stimulusSize=1000; % in pixels, 1000 = full-field stimulus
params.anesthOrAwake='anesth';
params.anesthType='iso';
params.LEDintensity=5;
params.peakWait=0.05; % Time in s 'til response onset/peak
fracSuppParams.baseline=[0.6 0.999];
% fracSuppParams.baseline=2+[0.6 0.999];
fracSuppParams.wait=0.04;

for i=1:length(orderedByDepthSpikes)
    orderedByDepthSpikes{i}=filtspikes(orderedByDepthSpikes{i},0,'fileInd',params.useFileInd);
end

if justExcludeAssigns==1 && doMU==0
    for i=1:length(criterionAssigns)
        currAssigns=criterionAssigns{i};
        currS=orderedByDepthSpikes{i};
        a=unique(currS.assigns);
        if isempty(criterionAssigns{i})
            newAssigns=a;
        else
            newAssigns=a(~ismember(a,currAssigns));
        end
        criterionAssigns{i}=newAssigns;
    end
elseif doMU==1
    for i=1:length(criterionAssigns)
        currAssigns=criterionAssigns{i};
        currS=orderedByDepthSpikes{i};
%         if ~isfield(currS,'assigns')
%             criterionAssigns{i}=nan;
%         end
        criterionAssigns{i}=nan;
    end
end

params.peakWindow=[params.stimulusOn(1)+params.peakWait params.stimulusOn(1)+(params.stimulusDuration/1000)-0.01];
fracSuppParams.ledWindow=[params.stimulusOn(1)+(params.stimulusDuration/1000)+fracSuppParams.wait params.stimulusOn(1)+(params.stimulusDuration/1000)+(params.LEDduration/1000)];

if spontExpt==0
    disp('DOING EVOKED');
else
    disp('DOING SPONT');
end

% Set up save directory for this mouse
saveDirName=[saveDir expt.name];
if exist(saveDirName,'dir')==0
    st=mkdir(saveDirName);
    if ~st
        disp('Could not create save directory');
    end
end

% Calculate information for criterion single units and save info. struct
assignsinfo.trode=[];
assignsinfo.original_assigns=[];
assignsinfo.event_channel=[];
assignsinfo.waveforms=[];
assignsinfo.waveformWidths=[];
criterionSpikes=[];
k=1;
for i=1:length(orderedByDepthSpikes)
    spikes=orderedByDepthSpikes{i};
    currTrode=Trode_numbers(Trode_orderByDepth(i));
    currAssigns=criterionAssigns{Trode_orderByDepth(i)};
    if isempty(currAssigns)
        continue
    end
    if isnan(currAssigns)
        if isempty(criterionSpikes)
            criterionSpikes=spikes;
        else
            criterionSpikes=concatSpikes_noAssigns(criterionSpikes,spikes);
        end
        continue
    end
    a=currAssigns;
    evChs=zeros(1,length(a));
    wwidths=zeros(1,length(a));
    for j=1:length(a)
        use_a_inds=spikes.assigns==a(j);
        theseEvChs=spikes.info.detect.event_channel(use_a_inds);
        topEvCh=mode(theseEvChs);
        evChs(j)=topEvCh;
%         disp(a(j));
%         if isnan(topEvCh)
%             disp('Unable to get spikes for unit');
%             disp(a(j));
%             continue
%         end
        wvfrmsEvCh=spikes.waveforms(:,:,topEvCh);
        wvfrm=mean(wvfrmsEvCh(use_a_inds,:),1);
        assignsinfo.waveforms(k,:)=wvfrm;
        k=k+1;
        wwidths(j)=classifyUnitWaveformWidth(wvfrm,3.2*10^-4,32000);
    end   
    currEvchs=orderUnitsbyDepth_getEvChs(spikes,a);
    currWidths=wwidths;
    if isempty(criterionSpikes)
        if dLGNdiscard==1
            criterionSpikes=filtspikes(spikes,0,'assigns',currAssigns(currWidths>=0.22*10^-3));
            disp('discard');
            disp(currAssigns(currWidths<0.22*10^-3));
        else
            criterionSpikes=filtspikes(spikes,0,'assigns',currAssigns);
        end
    else
        if dLGNdiscard==1
            criterionSpikes=concatSpikes_shiftAssigns(criterionSpikes,filtspikes(spikes,0,'assigns',currAssigns(currWidths>=0.22*10^-3)));
            disp('discard');
            disp(currAssigns(currWidths<0.22*10^-3));
        else
            criterionSpikes=concatSpikes_shiftAssigns(criterionSpikes,filtspikes(spikes,0,'assigns',currAssigns));
        end
    end
    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns)).*currTrode];
    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns];
    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs];
    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths];
end
if ~doMU==1
    [order,newassignsinfo]=orderUnitsByDepth(assignsinfo,Trode_numbers(Trode_orderByDepth));
    newassignsinfo.waveforms=assignsinfo.waveforms(order,:);
    newassignsinfo.waveformWidths=assignsinfo.waveformWidths(order);
    save([saveDirName '\criterionUnitsInfo.mat'],'newassignsinfo');
end
    
% Get silencing info for experiment and save, using all criterion units
% Also check window in which to calculate silencing
save([saveDirName '\paramsForSilencing.mat'],'params');
exptData=getSilencingInfoForExpt_passInParams(criterionSpikes,expt,params,checkLED,aliLED,dataDir);
save([saveDirName '\exptData.mat'],'exptData');

% Calculate fractional suppression for the summed unit PSTH and save
fracSupp=calcFractionalSupp_passInParams(exptData.xpoints,exptData.ypoints1,exptData.ypoints2,fracSuppParams,spontExpt);
save([saveDirName '\fracSupp_params.mat'],'fracSuppParams');
save([saveDirName '\fracSupp.mat'],'fracSupp');

if doMU==1
    return
end

% Calculate information for best single units and save info. struct
assignsinfo.trode=[];
assignsinfo.original_assigns=[];
assignsinfo.event_channel=[];
assignsinfo.waveforms=[];
assignsinfo.waveformWidths=[];
bestSpikes=[];
k=1;
for i=1:length(orderedByDepthSpikes)
    spikes=orderedByDepthSpikes{i};
    currTrode=Trode_numbers(Trode_orderByDepth(i));
    currAssigns=bestAssigns{Trode_orderByDepth(i)};
    if isempty(currAssigns)
        continue
    end
    a=currAssigns;
    evChs=zeros(1,length(a));
    wwidths=zeros(1,length(a));
    for j=1:length(a)
        use_a_inds=spikes.assigns==a(j);
        theseEvChs=spikes.info.detect.event_channel(use_a_inds);
        topEvCh=mode(theseEvChs);
        evChs(j)=topEvCh;
        wvfrmsEvCh=spikes.waveforms(:,:,topEvCh);
        wvfrm=mean(wvfrmsEvCh(use_a_inds,:),1);
        assignsinfo.waveforms(k,:)=wvfrm;
        k=k+1;
        wwidths(j)=classifyUnitWaveformWidth(wvfrm,3.2*10^-4,32000);
    end   
    currEvchs=orderUnitsbyDepth_getEvChs(spikes,a);
    currWidths=wwidths;
    if isempty(bestSpikes)
        bestSpikes=filtspikes(spikes,0,'assigns',currAssigns);
    else
        bestSpikes=concatSpikes_shiftAssigns(bestSpikes,filtspikes(spikes,0,'assigns',currAssigns));
    end
    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns)).*currTrode];
    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns];
    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs];
    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths];
end
[order,newassignsinfo]=orderUnitsByDepth(assignsinfo,Trode_numbers(Trode_orderByDepth));
newassignsinfo.waveforms=assignsinfo.waveforms(order,:);
newassignsinfo.waveformWidths=assignsinfo.waveformWidths(order);
save([saveDirName '\bestUnitsInfo.mat'],'newassignsinfo');

% Calculate silencing for individual units and save
[~,led1Spikes,led2Spikes,ns1,ns2,ass]=getUnitByUnitEffects(bestSpikes,spikes,unique(bestSpikes.assigns),params.useLEDcond{1},params.useLEDcond{2},params.useStimcond,params.useStimcond,fracSuppParams.ledWindow);
unitsFx.led1Spikes=led1Spikes;
unitsFx.led2Spikes=led2Spikes;
unitsFx.n1=ns1;
unitsFx.n2=ns2;
unitsFx.assigns=ass;
save([saveDirName '\unitsFx.mat'],'unitsFx');

% Save criterion units
save([saveDirName '\criterionSpikes.mat'],'criterionSpikes');

% Save best units
save([saveDirName '\bestSpikes.mat'],'bestSpikes');