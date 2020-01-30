function [psth,autocorr_out,powerOut,phaseOut]=measureAllUnitsResponseProperties(spikes,useAssigns,autocorr_window)

% Primarily for spontaneous activity in V1, where no importance of stimulus conditions

autocorr_out=[];
powerOut=[];
phaseOut=[];

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% useLED=[0 0.05 0.1515 0.2525 0.7575 5.00 5.05];
% useLED=[1.0000    1.0500    2.0000    2.0500    4.0000    4.0500    6.0000    6.0500    8.0000    8.0500   10.0000 ...
%         10.0500   12.0000   12.0500   14.0000   14.0500   16.0000   16.0500   18.0000   18.0500   20.0000   20.0500 ...
%         30.0000   30.0500   40.0000   40.0500   50.0000   50.0500   60.0000   60.0500];
% useLED=[freqs/100 freqs/100+5];
% useLED=sort(useLED);
% useLED=[0 0.1515 0.2525];
% useLED=[0 0.05 5.00 5.05];
% useLED=[0 5];
% useLED=[0 0.05 0.101 0.303 0.505 1.01 1.03 1.515 2.02 2.05 5 5.02 5 5.05 6 6.06 6.5 12.06 12.12 50.50];
% useLED=freqs+0.05;
useLED=sort([freqs freqs+0.05]);
% bin=1.2; % in ms
% bin=20; % in ms
bin=10; % in ms
% bin=1.5; % in ms
% trialDuration=5;
trialDuration=autocorr_window(2)-autocorr_window(1);
% useStimcond={[0:100000]};
% useStimcond={[1:8]};
useStimcond={[1.0000    1.0500    2.0000    2.0500    4.0000    4.0500    6.0000    6.0500    8.0000    8.0500   10.0000 ...
        10.0500   12.0000   12.0500   14.0000   14.0500   16.0000   16.0500   18.0000   18.0500   20.0000   20.0500 ...
        30.0000   30.0500   40.0000   40.0500   50.0000   50.0500   60.0000   60.0500]};

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
        useSpikes=prep_filtspikes(spikes,'stimcond',[],useStimcond{j});
        useSpikes=filtspikes(useSpikes,0,'assigns',allAssigns(i));
        %useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(single(currStimCon),single(useStimcond{j})) & ismember(single(currLEDCon),single(useLED)));
        else
            unitByUnitConsensus=currTrialCon(ismember(single(currStimCon),single(useStimcond{j})));
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