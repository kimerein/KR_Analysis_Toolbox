function [runningSU,bigFres]=bootstrapFxOfTemporalFreq(spikes,noThetaTrials,isNoTheta)

useLEDcond=[5.0100 ...
    5.0200    5.0400    5.0600    5.0800    5.1000    5.1200    5.1400    5.1600 ...
    5.1800    5.2000    5.3000    5.4000    5.5000    5.6000];
% Filt spikes by led condition and brain state
newass.new_assigns=unique(spikes.assigns);
tri=unique(spikes.trials);
useTrials=tri(noThetaTrials==isNoTheta);
spikes=filtspikes(spikes,0,'trials',useTrials);
whichTrialsPerLED=cell(1,length(useLEDcond));
for i=1:length(useLEDcond)
    spikes=makeTempField(spikes,'led',useLEDcond(i));
    whichTrialsPerLED{i}=unique(spikes.trials(spikes.temp==1));
    temp(i,:)=spikes.temp;
    temp1(i,:)=spikes.sweeps.temp;
    spikes.temp=[];
    spikes.sweeps.temp=[];
end
spikes.temp=sum(temp,1)>=1;
spikes.sweeps.temp=sum(temp1,1)>=1;
spikes=filtspikes(spikes,0,'temp',1);
useTrials=unique(spikes.trials);

% Bootstrap to get 95% confidence interval for all units
n=40;
nTrials=length(useTrials);
for i=1:n
    disp(i);
    % Sample with replacement for each LED condition
    currTrials=[];
    for j=1:length(useLEDcond)
        currTrials=[currTrials datasample(whichTrialsPerLED{j},length(whichTrialsPerLED{j}))];
    end
    [allSU,bigFres]=getFreqResponseInfo_forUnits_subFnctn([],[],makeSpikesStructFromSample(spikes,currTrials),newass,[],[]);
    if i==1
        runningSU=allSU;
    else
        for j=1:length(allSU)
            runningSU(j).p{length(runningSU(j).p)+1}=allSU(j).p{1};
            runningSU(j).responses{length(runningSU(j).responses)+1}=allSU(j).responses{1};
            runningSU(j).pStim{length(runningSU(j).pStim)+1}=allSU(j).pStim{1};
            runningSU(j).pSpont{length(runningSU(j).pSpont)+1}=allSU(j).pSpont{1};
        end
    end
end  

end

function newSpikes=makeSpikesStructFromSample(spikes,takeTrials)

takeTrials=sort(takeTrials);
newSpikes.led=[];
newSpikes.spiketimes=[];
newSpikes.trials=[];
newSpikes.assigns=[];
newSpikes.sweeps.fileInd=[];
newSpikes.sweeps.trials=[];
newSpikes.sweeps.led=[];

for i=1:length(takeTrials)
    currTrial=takeTrials(i); 
    subSpikes=filtspikes(spikes,0,'trials',currTrial);
    newSpikes.led=[newSpikes.led subSpikes.led];
    newSpikes.spiketimes=[newSpikes.spiketimes subSpikes.spiketimes];
    newSpikes.trials=[newSpikes.trials i*ones(size(subSpikes.trials))];
    newSpikes.assigns=[newSpikes.assigns subSpikes.assigns];
    newSpikes.sweeps.fileInd=[newSpikes.sweeps.fileInd unique(subSpikes.sweeps.fileInd(ismember(subSpikes.sweeps.trials,currTrial)))];
    newSpikes.sweeps.trials=[newSpikes.sweeps.trials i];
    newSpikes.sweeps.led=[newSpikes.sweeps.led unique(subSpikes.sweeps.led(ismember(subSpikes.sweeps.trials,currTrial)))];
end

end

function [allSU,bigFres]=getFreqResponseInfo_forUnits_subFnctn(expt,orderedByDepthSpikes,alreadySpikes,newass,allSpecs,f_forAllSpecs)
% Start with non-filtered spikes, so can read event_channel out of struct
% without having to recalculate

justExcludeAssigns=0;
doMU=0;
filtBySpecs=0;

% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked dLGN stim. duration 0 ms\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont dLGN\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Evoked LP\';
% saveDir='W:\Analysis Computer\Thalamus Silencing Across Mice\Spont LP\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Evoked Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Spont Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Static Silencing Across Mice\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\Chronux SU Thal No Trial-Average\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\Chronux Cortex SU SS\';
% saveDir='W:\Analysis Computer\Across Mice Finals\Frequency Response\Chronux Cx SU Habituation\';
% dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';

saveDir='F:\Figures\Different Temporal Frequencies\M384\newer theta thresh\theta LED\';
dataDir='W:\New Acquisition Computer\';

% Fs=32000;
Fs=25000;
spontExpt=0; % 0 if evoked silencing experiment
dLGNdiscard=0; % don't use units with half-width <0.22 ms in frac. supp. and PSTH

Trode_numbers=[5 6 9]; % Must be in order of depth
Trode_orderByDepth=[1 2 3];

criterionAssigns{1}=[15]; % This matches Trode_number
criterionAssigns{2}=[23 83 94 97 99];
criterionAssigns{3}=[17 35 38 39 48 49 50 60 67 71 72 78 82 96 106 124];
criterionAssigns{4}=[77 89 97 101 115];
% criterionAssigns{5}=[10 11 15 22 24 32 34 39 41 42 44 46 53 57 58 59 64 67 68 69 70 71 72 78 80 82 83 87 91 93 98];

% bestAssigns{1}=[]; % This matches Trode_number
% bestAssigns{2}=[];
% bestAssigns{3}=[];
% bestAssigns{4}=[];

params.freqBand=[4 8]; % in Hz
params.doSingleTrials=0;
doSingleTrials=params.doSingleTrials;
params.powerThresh=0.1;
params.useFileInd=[1:100000];
params.useStimcond=[1:128];
params.useTrigger=1:12;
% params.useLEDcond={[1.0000 2.0000 4.0000 6.0000 8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000 30.0000 40.0000 50.0000 60.0000]};
% params.useLEDcond={[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050]};
% params.useLEDcond={[1.010 2.010 4.010 6.010 8.010 10.010 12.010 14.010 16.010 18.010 20.010 30.010 40.010 50.010 60.010 ...
%     1.030 2.030 4.030 6.030 8.030 10.030 12.030 14.030 16.030 18.030 20.030 30.030 40.030 50.030 60.030]};
% params.useLEDcond={[0.0100    0.0200    0.0400    0.0600    0.0800    0.1000    0.1200    0.1400 ...
%     0.1600    0.1800    0.2000    0.3000    0.4000    0.5000    0.6000]};   
params.useLEDcond={[5.0100 ...
    5.0200    5.0400    5.0600    5.0800    5.1000    5.1200    5.1400    5.1600 ...
    5.1800    5.2000    5.3000    5.4000    5.5000    5.6000]};
% params.useLEDcond={[1.0000 2.0000 4.0000 6.0000 8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000 30.0000 40.0000 50.0000 60.0000],[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050]};
% params.useLEDcond={[1 1.050 2 2.050 4 4.050 6 6.050 8 8.050 10 10.050 12 12.050 14 14.050 16 16.050 18 18.050 20 20.050 30 30.050 40 40.050 50 50.050 60 60.050]};
params.trialDuration=4; % in s
params.stimulusOn=[1 3]; % in s relative to trial onset
params.stimulusType='Sine';
params.LEDonset=500; % ms into trial for LED onset
params.LEDduration=3000; % in ms, duration of LED pulse -- maybe short to exclude photo-artifact
params.baseWindow=[0 1];
params.analysisType='cross-corr'; 
params.anesthOrAwake='awake';
params.anesthType='na';
params.LEDintensity=5;


% Set up save directory for this mouse
% saveDirName=[saveDir expt.name];
% if exist(saveDirName,'dir')==0
%     st=mkdir(saveDirName);
%     if ~st
%         disp('Could not create save directory');
%     end
% end
% save([saveDirName '\params.mat'],'params');

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
    assignsinfo.new_assigns=[];
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
                wwidths(j)=classifyUnitWaveformWidth(wvfrm,3.2*10^-4,32000);
                if dLGNdiscard==1
                    if wwidths(j)>=0.22*10^-3
                        assignsinfo.waveforms(k,:)=wvfrm;
                        k=k+1;
                    end
                else
                    assignsinfo.waveforms(k,:)=wvfrm;
                    k=k+1;
                end
            end
            currEvchs=orderUnitsbyDepth_getEvChs(spikes,a);
            currWidths=wwidths;
            if isempty(criterionSpikes)
                if dLGNdiscard==1
                    criterionSpikes=filtspikes(spikes,0,'assigns',currAssigns(currWidths>=0.22*10^-3));
                    disp('discard');
                    disp(currAssigns(currWidths<0.22*10^-3));
                    assignsinfo.new_assigns=currAssigns(currWidths>=0.22*10^-3);
                    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns(currWidths>=0.22*10^-3))).*currTrode];
                    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns(currWidths>=0.22*10^-3)];
                    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs(currWidths>=0.22*10^-3,:)];
                    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths(currWidths>=0.22*10^-3)];
                else
                    criterionSpikes=filtspikes(spikes,0,'assigns',currAssigns);
                    assignsinfo.new_assigns=currAssigns;
                    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns)).*currTrode];
                    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns];
                    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs];
                    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths];
                end
            else
                if dLGNdiscard==1
                    criterionSpikes=concatSpikes_shiftAssigns(criterionSpikes,filtspikes(spikes,0,'assigns',currAssigns(currWidths>=0.22*10^-3)));
                    disp('discard');
                    disp(currAssigns(currWidths<0.22*10^-3));
                    assignsinfo.new_assigns=[assignsinfo.new_assigns max(assignsinfo.new_assigns)+currAssigns(currWidths>=0.22*10^-3)];
                    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns(currWidths>=0.22*10^-3))).*currTrode];
                    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns(currWidths>=0.22*10^-3)];
                    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs(currWidths>=0.22*10^-3,:)];
                    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths(currWidths>=0.22*10^-3)];
                else
                    criterionSpikes=concatSpikes_shiftAssigns(criterionSpikes,filtspikes(spikes,0,'assigns',currAssigns));
                    assignsinfo.new_assigns=[assignsinfo.new_assigns max(assignsinfo.new_assigns)+currAssigns];
                    assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns)).*currTrode];
                    assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns];
                    assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs];
                    assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths];
                end
            end 
        end
        [order,newassignsinfo]=orderUnitsByDepth(assignsinfo,Trode_numbers(Trode_orderByDepth));
        newassignsinfo.waveforms=assignsinfo.waveforms(order,:);
        newassignsinfo.waveformWidths=assignsinfo.waveformWidths(order);
        newassignsinfo.new_assigns=assignsinfo.new_assigns(order);
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
%     criterionSpikes=filtspikes(alreadySpikes,0,'fileInd',params.useFileInd);
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
      
% Get SU (criterionSpikes) frequency response
if ~isempty(newass)
    a=newass.new_assigns;
else
    a=newassignsinfo.new_assigns;
end
if doSingleTrials==0
    for k=1:length(a)
%         disp(k);
        for i=1:length(params.useLEDcond)
            ledValue=params.useLEDcond{i};
%             temp=[];
%             temp1=[];
%             for j=1:length(ledValue)
%                 criterionSpikes=makeTempField(criterionSpikes,'led',ledValue(j));
%                 temp(j,:)=criterionSpikes.temp;
%                 temp1(j,:)=criterionSpikes.sweeps.temp;
%                 criterionSpikes.temp=[];
%                 criterionSpikes.sweeps.temp=[];
%             end
%             criterionSpikes.temp=sum(temp,1)>=1;
%             criterionSpikes.sweeps.temp=sum(temp1,1)>=1;
%             currSpikes=filtspikes(criterionSpikes,0,'temp',1);
            currSpikes=criterionSpikes;
            [freqs,p,responses,pStim,pSpont,bigFres]=getUnitFreqResponse_matrix(currSpikes,a(k),ledValue);
            SUcrosscorr.freqs{i}=freqs;
            SUcrosscorr.p{i}=p;
            SUcrosscorr.responses{i}=responses;
            SUcrosscorr.pStim{i}=pStim;
            SUcrosscorr.pSpont{i}=pSpont;
            SUcrosscorr.assigns{i}=a(k);
        end
        allSU(k)=SUcrosscorr;
    end
%     save([saveDirName '\allSU' num2str(i) '.mat'],'allSU');
%     save([saveDirName '\bigFres' num2str(i) '.mat'],'bigFres');
else
    for k=1:length(a)
        disp(k);
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
            [freqs,p,responses,pStim,pSpont,bigFres]=getUnitFreqResponse_matrix_singleTrials(currSpikes,a(k),ledValue);
            SUcrosscorr.freqs{i}=freqs;
            SUcrosscorr.p{i}=p;
            SUcrosscorr.responses{i}=responses;
            SUcrosscorr.pStim{i}=pStim;
            SUcrosscorr.pSpont{i}=pSpont;
            SUcrosscorr.assigns{i}=a(k);
        end
        allSU(k)=SUcrosscorr;
    end
%     save([saveDirName '\allSUsingleTrials' num2str(i) '.mat'],'allSU');
%     save([saveDirName '\bigFres' num2str(i) '.mat'],'bigFres');
end

% Calculate information for best single units and save info. struct
% if doMU==0
%     assignsinfo.trode=[];
%     assignsinfo.original_assigns=[];
%     assignsinfo.event_channel=[];
%     assignsinfo.waveforms=[];
%     assignsinfo.waveformWidths=[];
%     bestSpikes=[];
%     k=1;
%     for i=1:length(orderedByDepthSpikes)
%         spikes=orderedByDepthSpikes{i};
%         currTrode=Trode_numbers(Trode_orderByDepth(i));
%         currAssigns=bestAssigns{Trode_orderByDepth(i)};
%         if isempty(currAssigns)
%             continue
%         end
%         a=currAssigns;
%         evChs=zeros(1,length(a));
%         wwidths=zeros(1,length(a));
%         for j=1:length(a)
%             use_a_inds=spikes.assigns==a(j);
%             theseEvChs=spikes.info.detect.event_channel(use_a_inds);
%             topEvCh=mode(theseEvChs);
%             evChs(j)=topEvCh;
%             wvfrmsEvCh=spikes.waveforms(:,:,topEvCh);
%             wvfrm=mean(wvfrmsEvCh(use_a_inds,:),1);
%             assignsinfo.waveforms(k,:)=wvfrm;
%             k=k+1;
%             wwidths(j)=classifyUnitWaveformWidth(wvfrm,3.2*10^-4,32000);
%         end
%         currEvchs=orderUnitsbyDepth_getEvChs(spikes,a);
%         currWidths=wwidths;
%         if isempty(bestSpikes)
%             bestSpikes=filtspikes(spikes,0,'assigns',currAssigns);
%         else
%             bestSpikes=concatSpikes_shiftAssigns(bestSpikes,filtspikes(spikes,0,'assigns',currAssigns));
%         end
%         assignsinfo.trode=[assignsinfo.trode ones(1,length(currAssigns)).*currTrode];
%         assignsinfo.original_assigns=[assignsinfo.original_assigns currAssigns];
%         assignsinfo.event_channel=[assignsinfo.event_channel; currEvchs];
%         assignsinfo.waveformWidths=[assignsinfo.waveformWidths currWidths];
%     end
%     [order,newassignsinfo]=orderUnitsByDepth(assignsinfo,Trode_numbers(Trode_orderByDepth));
%     newassignsinfo.waveforms=assignsinfo.waveforms(order,:);
%     newassignsinfo.waveformWidths=assignsinfo.waveformWidths(order);
%     save([saveDirName '\bestUnitsInfo.mat'],'newassignsinfo');
%     
%     % Calculate silencing for individual units and save
%     for i=1:length(params.useLEDcond)
% %         currLED=params.useLEDcond{i};
% %         currSpikes=filtspikes(bestSpikes,0,'led',currLED);
% %         [freqs,p,responses,pStim,pSpont]=getAllFreqResponse_matrix(bestSpikes);
% %         SUcrosscorr.freqs=freqs;
% %         SUcrosscorr.p=p;
% %         SUcrosscorr.responses=responses;
% %         SUcrosscorr.pStim=pStim;
% %         SUcrosscorr.pSpont=pSpont;
% %         save([saveDirName '\SUcrosscorr_' num2str(i) '.mat'],'SUcrosscorr');
%     end
% end

% Save criterion units
% save([saveDirName '\criterionSpikes.mat'],'criterionSpikes');

% Save best units
% if doMU==0
%     save([saveDirName '\bestSpikes.mat'],'bestSpikes');
% end

end