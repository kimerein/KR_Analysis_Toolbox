function [output,psth]=find_putative_ArchT_units(spikes,psth,LEDbySweep,ledConds,taketri,ledVal,trialDuration)

% direct_window=35; % in ms, the time from LED onset below which a 
% significant suppression indicates ArchT expression
% baseline_window=[0 3]; % in seconds from trial onset
% stimWindow=[4 6.5]; % in seconds from trial onset

%  For V1 units
direct_window=6000; % in ms, the time from LED onset below which a 
% significant suppression indicates ArchT expression
baseline_window=[0 3]; % in seconds from trial onset
stimWindow=[4 6.5]; % in seconds from trial onset

% Analyze the following:
% 1. Amplitude of the change from baseline during the first (default) 35-45 ms after
% LED onset
% 2. P-val of the change from baseline during this same time window
% 3. Suppression of the evoked response
% 4. Waveform half-width-at-half-max

if isempty(psth)
    spikes=filtspikes(spikes,0,'trials',taketri);
    spikes.trials=spikes.trials-min(spikes.trials)+1;
    spikes.sweeps.trials=spikes.sweeps.trials-min(spikes.sweeps.trials)+1;
    [psth]=measureAllUnitsResponseProperties(spikes,unique(spikes.assigns),[0 14.5]);
end

if ~isempty(LEDbySweep)
    if iscell(LEDbySweep)
        LEDbySweep=LEDbySweep{1};
    end
    
    if size(LEDbySweep,1)~=size(psth.psths{1},1)
        error('Sizes of LEDbySweep and psth do not match');
    end
    if ~isequal(ledConds,psth.unitLED{1})
        error('LED conditions for LEDbySweep and psth do not match');
    end
    
    % Align psth to led
    times=linspace(0,trialDuration,size(LEDbySweep,2));
    [psth,ledConds,ledStart]=alignPSTHtoLED(psth,LEDbySweep,ledConds,ledVal,times,[]);
else
    ledStart=3.2;
end
times=psth.t;

% 1. Amplitude of the change from baseline during the first (default) 35-45 ms after
% LED onset
direct_window=direct_window/1000; % convert to seconds
amp_change=nan(1,length(psth.psths));
% 2. P-val of the change from baseline during this same time window
pvals=nan(1,length(psth.psths));
% 3. Suppression of the evoked response
visev_change=nan(1,length(psth.psths)); % amplitude change during visual stimulus
visev_pvals=nan(1,length(psth.psths)); % p-val of change during visual stimulus
for i=1:length(psth.psths)
    currp=psth.psths{i};
    l=psth.unitLED{1};
    bases=currp(ismember(single(l),single(ledVal)),times>=baseline_window(1) & times<baseline_window(2));
    during_led=currp(ismember(single(l),single(ledVal)),times>=ledStart & times<ledStart+direct_window);
    spont_no_led=currp(~ismember(single(l),single(ledVal)),times>=ledStart & times<ledStart+direct_window);
    evoked_noLED=nanmean(currp(~ismember(single(l),single(ledVal)),times>=stimWindow(1) & times<stimWindow(2)),2);
    evoked_LED=nanmean(currp(ismember(single(l),single(ledVal)),times>=stimWindow(1) & times<stimWindow(2)),2);
%     amp_change(i)=nanmean(nanmean(during_led,2)-nanmean(bases,2));
    amp_change(i)=nanmean(nanmean(spont_no_led,2))-nanmean(nanmean(during_led,2));
    pvals(i)=ranksum(nanmean(spont_no_led,2),nanmean(during_led,2));
    visev_change(i)=nanmean(evoked_noLED)-nanmean(evoked_LED);
    visev_pvals(i)=ranksum(evoked_noLED,evoked_LED);
end

cmap=colormap(jet(20));
color_edges=fliplr(0:0.05:0.95);
color_ind=0;

figure();
output.amp_change=amp_change;
output.pvals=pvals;
output.visev_change=visev_change;
output.colors=nan(length(pvals),3);
output.visev_pvals=visev_pvals;
for i=1:length(pvals)
    p=pvals(i);
    color_ind=find(p>color_edges,1,'first');
    if isempty(color_ind)
        color_ind=1;
    end
%     scatter(amp_change(i),visev_change(i),[],cmap(color_ind,:));
    scatter(1-p,visev_change(i),[],cmap(color_ind,:));
    output.colors(i,:)=cmap(color_ind,:);
    hold on;
end
colorbar;
xlabel('amp change of suppression at led onset');
ylabel('change in activity during visual stimulus');
title('Want positive x and y');

end

function [psth,autocorr_out,powerOut,phaseOut]=measureAllUnitsResponseProperties(spikes,useAssigns,autocorr_window)

% Primarily for spontaneous activity in V1, where no importance of stimulus conditions

autocorr_out=[];
powerOut=[];
phaseOut=[];

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
useLED=[0 0.05 0.1515 0.2525 0.7575 5.00 5.05];
% useLED=[0 0.1515 0.2525];
% useLED=[0 0.05 5.00 5.05];
% useLED=freqs+0.05;
% useLED=sort([freqs freqs+0.05]);
% bin=1.2; % in ms
% bin=20; % in ms
bin=5; % in ms
% trialDuration=5;
trialDuration=autocorr_window(2)-autocorr_window(1);
useStimcond={[1:100000]};

% Get trial by trial PSTHs for units
[x,psths_t,unitTrials,unitStimcond,unitLED]=getTrialByTrialUnitPSTH_sub(spikes,useAssigns,useLED,bin,trialDuration,useStimcond);
x=[x max(x)+(x(2)-x(1))];
psth.t=x;
psth.psths=psths_t;
psth.unitTrials=unitTrials;
psth.unitStimcond=unitStimcond;
psth.unitLED=unitLED;

% Get auto-correlation of each unit, each trial
autocorrs=cell(1,length(psths_t));
autocorr_lags=[];
% autocorr_window=[1.2 2.7];
% autocorr_window=[3 5];
% autocorr_window=[3 3.5];
% autocorr_window=[1.5 3];
% autocorr_window=[0.25 1.75];
for i=1:length(psths_t)
    disp('Cell #');
    disp(i);
    curr_psth=psths_t{i};
    % Get auto-correlation of each trial
    for j=1:size(curr_psth,1)
        c=curr_psth(j,:);
%         takec=c((x>=0 & x<=0.25) | (x>=1.75 & x<=2));
%         [r,lags]=xcorr(takec,takec,floor(2/(bin/1000)),'coeff');
        [r,lags]=xcorr(c(x>=autocorr_window(1) & x<=autocorr_window(2)),c(x>=autocorr_window(1) & x<=autocorr_window(2)),floor(2/(bin/1000)),'coeff');
        if j==1
            autocorr_lags=lags;
            temp_autocorr=zeros(size(curr_psth,1),length(r));
        end
        if any(isnan(r))
            r=zeros(size(r));
        end
        temp_autocorr(j,:)=r;
    end
    autocorrs{i}=temp_autocorr;
end
autocorr_out.autocorrs=autocorrs;
autocorr_out.lags=autocorr_lags;

return

% Get each unit's frequency, phase and amplitude as a function of time and
% LED condition
params.tapers=[5 9];
params.Fs=1/(bin/1000);
params.fpass=[1 100];
params.trialave=0;
% Assumption 1: Look for persistent cell assemblies of more than 500 ms
movingwin=[0.5 0.01];
doFreqsForPhase=7;
paramsForPhase=params;

powerSpectra=cell(1,length(psths_t));
phaseSpectra=cell(1,length(psths_t));
for i=1:length(psths_t)
    disp('Cell #');
    disp(i);
    curr_psth=psths_t{i};
    % Power as a function of time and frequency
    [S,t,f]=mtspecgrampb(curr_psth',movingwin,params);
    if i==1
        powerSpectra_times=t;
        powerSpectra_freqs=f;
    end
    powerSpectra{i}=S; % time x freq x trials
    % Phase as a function of time and frequency
    for j=1:length(doFreqsForPhase)
        disp('phase analysis');
        disp(j);
        currFreq=doFreqsForPhase(j);
        referentSine=repmat(sin(2*pi*currFreq.*x),size(curr_psth,1),1);
        paramsForPhase.fpass=[currFreq-0.5 currFreq+0.5];
        [C,phi,~,~,~,t,f]=cohgramcpb(referentSine',curr_psth',movingwin,params);
        Cav=nanmean(C,2);
        if j==1
            currPhaseSpectra=nan(size(C,1),length(doFreqsForPhase),size(C,3));
            phaseSpectra_times=t;
            phaseSpectra_freqs=f;
        end
        currPhaseSpectra(:,j,:)=Cav;
    end
    phaseSpectra{i}=currPhaseSpectra; % time x freq x trials
end
powerOut.powerSpectra=powerSpectra;
powerOut.t=powerSpectra_times;
powerOut.f=powerSpectra_freqs;
phaseOut.phaseSpectra=phaseSpectra;
phaseOut.t=phaseSpectra_times;
phaseOut.f=phaseSpectra_freqs;

end

function [x,psths_t,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED]=getTrialByTrialUnitPSTH_sub(spikes,allAssigns,useLED,bin,trialDuration,useStimcond)

% Get trials for each unit
unitByUnitTrials=cell(1,length(allAssigns));
unitByUnitStimcond=cell(1,length(allAssigns));
unitByUnitLED=cell(1,length(allAssigns));
tt=unique(spikes.sweeps.trials);
if tt(1)>1 && length(spikes.sweeps.trials)==length(spikes.sweeps.stimcond)
    spikes.sweeps.trials=spikes.sweeps.trials-tt(1)+1;
    spikes.trials=spikes.trials-tt(1)+1;
end
% if any(isnan(spikes.sweeps.led))
%     spikes.sweeps.led(isnan(spikes.sweeps.led))=-10;
%     spikes.led(isnan(spikes.led))=-10;
% end
% if any(isnan(spikes.sweeps.stimcond))
%     spikes.sweeps.stimcond(isnan(spikes.sweeps.stimcond))=-10;
%     spikes.stimcond(isnan(spikes.stimcond))=-10;
% end

% Changed to deal with nans 7/22/15
for i=1:length(allAssigns)
    strials=spikes.sweeps.trials(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    sstim=spikes.sweeps.stimcond(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    sled=spikes.sweeps.led(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    [unitByUnitTrials{i},indsin]=unique(strials);
    unitByUnitStimcond{i}=sstim(indsin);
    unitByUnitLED{i}=sled(indsin); 
end

% for i=1:length(allAssigns)
%     unitByUnitTrials{i}=unique(spikes.sweeps.trials);
%     unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
%     unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
% end

psths_t=cell(length(allAssigns),length(useStimcond));
x=[];
for i=1:length(allAssigns)
    for j=1:length(useStimcond)
        useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(single(currLEDCon),single(useLED)));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        if length(unitByUnitConsensus)~=length(unitByUnitTrials{i})
            disp('stop here');
        end
        [~,~,~,x,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        if size(psths_t{i,j},1)~=length(unitByUnitTrials{i})
            disp('stop');
        end
    end
end

end


    
    



