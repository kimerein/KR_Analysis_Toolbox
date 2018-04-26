function getFreqResponseInfo(expt,orderedByDepthSpikes,alreadySpikes,allSpecs,f_forAllSpecs)
% Start with non-filtered spikes, so can read event_channel out of struct
% without having to recalculate

justExcludeAssigns=0;
doMU=1;
filtBySpecs=0;

% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN stim. duration 0 ms\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked LP\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont LP\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Evoked Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Spont Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Static Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\Summed Single Units Chronux\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\SS Cortex MU Chronux Anesth\';
saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\Chronux Cx SU Habituation\';
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
Fs=32000;
spontExpt=0; % 0 if evoked silencing experiment
dLGNdiscard=0; % don't use units with half-width <0.22 ms in frac. supp. and PSTH

Trode_numbers=[1]; % Must be in order of depth
Trode_orderByDepth=[1];

% criterionAssigns{1}=[]; % This matches Trode_number
% criterionAssigns{2}=[];
% criterionAssigns{3}=[];
% criterionAssigns{4}=[];

% bestAssigns{1}=[]; % This matches Trode_number
% bestAssigns{2}=[];
% bestAssigns{3}=[];
% bestAssigns{4}=[];

params.freqBand=[4 8]; % in Hz
params.powerThresh=0.1;
params.useFileInd=[1:10000];
params.useStimcond=[1:128];
params.useTrigger=1:12;
% params.useLEDcond={[1.0000 2.0000 4.0000 6.0000 8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000 30.0000 40.0000 50.0000 60.0000]};
params.useLEDcond={[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050]};
params.trialDuration=4; % in s
params.stimulusOn=[1.4 3]; % in s relative to trial onset
params.stimulusType='Sine';
params.LEDonset=500; % ms into trial for LED onset
params.LEDduration=3000; % in ms, duration of LED pulse -- maybe short to exclude photo-artifact
params.baseWindow=[0 1];
params.analysisType='auto-corr'; 
params.anesthOrAwake='anesth';
params.anesthType='iso';
params.LEDintensity=0;
% Set up save directory for this mouse
saveDirName=[saveDir expt.name];
if exist(saveDirName,'dir')==0
    st=mkdir(saveDirName);
    if ~st
        disp('Could not create save directory');
    end
end
save([saveDirName '\params.mat'],'params');

if isempty(alreadySpikes)
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
    end
    
    if spontExpt==0
        disp('DOING EVOKED');
    else
        disp('DOING SPONT');
    end
    
    % Calculate information for criterion single units and save info. struct
    assignsinfo.trode=[];
    assignsinfo.original_assigns=[];
    assignsinfo.event_channel=[];
    assignsinfo.waveforms=[];
    assignsinfo.waveformWidths=[];
    criterionSpikes=[];
    k=1;
    if doMU==0
        for i=1:length(orderedByDepthSpikes)
            spikes=orderedByDepthSpikes{i};
            currTrode=Trode_numbers(Trode_orderByDepth(i));
            currAssigns=criterionAssigns{Trode_orderByDepth(i)};
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
        [order,newassignsinfo]=orderUnitsByDepth(assignsinfo,Trode_numbers(Trode_orderByDepth));
        newassignsinfo.waveforms=assignsinfo.waveforms(order,:);
        newassignsinfo.waveformWidths=assignsinfo.waveformWidths(order);
        save([saveDirName '\criterionUnitsInfo.mat'],'newassignsinfo');
    else
        for i=1:length(orderedByDepthSpikes)
            spikes=orderedByDepthSpikes{i};
            if isempty(criterionSpikes)
                criterionSpikes=spikes;
            else
                criterionSpikes=concatSpikes_noAssigns(criterionSpikes,spikes);
            end
        end
    end
else
    criterionSpikes=alreadySpikes;
end

if filtBySpecs==1
    if length(unique(criterionSpikes.trials))~=size(allSpecs,1)
        disp('mismatch in number of trials in allSpecs');
    else
        useTheseTrials=max(allSpecs(:,f_forAllSpecs>=params.freqBand(1) & f_forAllSpecs<=params.freqBand(2)),[],2)<params.powerThresh;
        ttt=unique(criterionSpikes.trials);
        criterionSpikes=filtspikes(criterionSpikes,0,'trials',ttt(useTheseTrials));
    end
end
      
% Get MU (criterionSpikes) frequency response
for i=1:length(params.useLEDcond)
    ledValue=params.useLEDcond{i};
    temp=[];
    temp1=[];
    for j=1:length(ledValue)
        criterionSpikes=makeTempField(criterionSpikes,'led',ledValue(j));
        temp(j,:)=criterionSpikes.temp;
        temp1(j,:)=criterionSpikes.sweeps.temp;
        criterionSpikes.temp=[];
        criterionSpikes.sweeps.temp=[];
    end
    criterionSpikes.temp=sum(temp,1)>=1;
    criterionSpikes.sweeps.temp=sum(temp1,1)>=1;
    currSpikes=filtspikes(criterionSpikes,0,'temp',1);
    [freqs,p,responses,pStim,pSpont,bigFres]=getUnitFreqResponse_matrix(currSpikes,[],ledValue);
    MUcrosscorr.freqs{i}=freqs;
    MUcrosscorr.p{i}=p;
    MUcrosscorr.responses{i}=responses;
    MUcrosscorr.pStim{i}=pStim;
    MUcrosscorr.pSpont{i}=pSpont;
end
save([saveDirName '\MUcrosscorr_' num2str(i) '.mat'],'MUcrosscorr');
save([saveDirName '\responseFrequencies.mat'],'bigFres');

% Calculate information for best single units and save info. struct
if doMU==0
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
    for i=1:length(params.useLEDcond)
%         currLED=params.useLEDcond{i};
%         currSpikes=filtspikes(bestSpikes,0,'led',currLED);
%         [freqs,p,responses,pStim,pSpont]=getAllFreqResponse_matrix(bestSpikes);
%         SUcrosscorr.freqs=freqs;
%         SUcrosscorr.p=p;
%         SUcrosscorr.responses=responses;
%         SUcrosscorr.pStim=pStim;
%         SUcrosscorr.pSpont=pSpont;
%         save([saveDirName '\SUcrosscorr_' num2str(i) '.mat'],'SUcrosscorr');
    end
end

% Save criterion units
% save([saveDirName '\criterionSpikes.mat'],'criterionSpikes');

% Save best units
if doMU==0
    save([saveDirName '\bestSpikes.mat'],'bestSpikes');
end