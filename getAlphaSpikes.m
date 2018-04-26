function [spikes,LFPbySweep]=getAlphaSpikes(expt,spikes,LFPbySweep,alphaPower)
% function [spikes,LFPbySweep]=getAlphaSpikes(expt,spikes,LFPbySweep,alphaPower)
% If LFPbySweep is empty, will read in and save LFP data from raw files
% If LFPbySweep is passed in, will use this value
% If alphaPower is empty, will re-make spectrograms
% If alphaPower is passed in, will use this value

%% For Lea
% Define analysis parameters here!
fileInd=[70:148]; % files to include in analysis
alphaBand=[4 8]; % frequency band considered Alpha, in Hz
% alphaBand=[4 10]; % frequency band considered Alpha, in Hz
saveDir='W:\Analysis Computer\Getting Alpha Spikes\'; % output directory for saving data
% dataDir='W:\New Acquisition Computer\'; % directory containing raw data
dataDir='W:\New Acquisition Computer\';
useChannelNumber=3; % use this channel from top of probe for LFP
trialDuration=3.5; % duration of a single trial in seconds
allStimcond=1:12; % which stimconds to include in analysis
exampleTrial=800;
Fs=32000;
% Fs=1000;
downSampFactor=10;
% downSampFactor=1;

% Change alpha thresh to pick out alpha states
% alphaThresh=40000; 
alphaThresh=60;
upperBand=[30 100];
lowerBand=alphaBand;
smoothBin=1500; % use if LFP
% smoothBin=100; % use if spiking

%%

if exampleTrial>length(fileInd)*12
    disp('Choose an example trial within fileInd');
    return
end

spikes=filtspikes(spikes,0,'fileInd',fileInd);

% Step 1: Read in LFP and find times with high alpha power
if ~isempty(LFPbySweep)
    if isempty(alphaPower)
        [UPs,alphaPower,ledToReturn,LFPbySweep]=find_alpha_wrapper(expt,useChannelNumber,fileInd,spikes,LFPbySweep,dataDir,alphaThresh,exampleTrial,Fs,downSampFactor,lowerBand,upperBand,smoothBin);
        save([saveDir 'alphaPower.mat'],'alphaPower');
    else
        ledToReturn=[];
        UPs=[];
    end
else
    [UPs,alphaPower,ledToReturn,LFPbySweep]=find_alpha_wrapper(expt,useChannelNumber,fileInd,spikes,[],dataDir,alphaThresh,exampleTrial,Fs,downSampFactor,lowerBand,upperBand,smoothBin);
    save([saveDir 'LFPbySweep.mat'],'LFPbySweep');
    save([saveDir 'alphaPower.mat'],'alphaPower');
end

% Step 2: Get alpha power-defined alpha states
x=linspace(0,trialDuration,length(alphaPower{1}));
alphaStates=returnAlpha(x,alphaPower,alphaThresh);
save([saveDir 'alphaStates.mat'],'alphaStates');

% Step 3: Determine whether spikes occur during alpha
if length(alphaStates)~=length(unique(spikes.sweeps.trials))
    disp('Trial numbers do not match');
    spikes=[];
else
    spikes.alpha=zeros(size(spikes.trials));
    spikeTrials=unique(spikes.sweeps.trials);
    for i=1:length(spikeTrials)
        currTrial=spikeTrials(i);
        currStates=alphaStates{i};
        checkSpiketimes=spikes.spiketimes(spikes.trials==currTrial);
        useInds=find(spikes.trials==currTrial);
        for j=1:length(checkSpiketimes)
            currSpiketime=checkSpiketimes(j);
            inAlphaState=0;
            for k=1:size(currStates,1)
                if currSpiketime>=currStates(k,1) & currSpiketime<=currStates(k,2)
                    inAlphaState=1;
                    break
                end
            end
            spikes.alpha(useInds(j))=inAlphaState;
        end
    end
end
    

end

function [UPstates,powerRatio,ledToReturn,LFPbySweep]=find_alpha_wrapper(expt,useChannelNumber,fileInd,spikes,LFPbySweep,dataDir,alphaThresh,exampleTrial,Fs,downSampFactor,lowerBand,upperBand,smoothBin)

powerRatio=[];
ledToReturn=[];
UPstates=[];

% lowerBand=[4 8];
% upperBand=[30 100];

% Fs=32000;
% downSampFactor=10;
noLEDcond=nan;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
usePhysCh_ind=useChannelNumber;

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames); 

% Get LFP
if isempty(LFPbySweep)
    LFPbySweep=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(usePhysCh_ind));
end
Fs=Fs/downSampFactor;
if any(size(LFPbySweep{1})==0)
    disp('No data in this daq file.');
end

% Match acquired LFP sweeps to LED conditions
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    return
end

% Get stimconds
stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));

temp=LFPbySweep{1};
subLFPbySweep=LFPbySweep;
if ~isnan(noLEDcond)
    subLFPbySweep=temp(ismember(ledConds,noLEDcond),:);
    subSpikes=filtspikes(spikes,0,'led',noLEDcond);
    ledToReturn=ledConds(ismember(ledConds,noLEDcond));
else
    subLFPbySweep=temp;
    subSpikes=spikes;
    ledToReturn=nan;
end
disp('These numbers should be the same');
disp(size(subLFPbySweep,1));
disp(length(unique(subSpikes.sweeps.trials)));

[UPstates,powerRatio]=find_alpha_states(subLFPbySweep,Fs,lowerBand,upperBand,alphaThresh,exampleTrial,smoothBin);
end

function [UPstates,powerRatio]=find_alpha_states(LFPbySweep,LFP_Fs,lowerBand,upperBand,alphaThresh,exampleTrial,smoothBin)

global dontShowSpecgram
dontShowSpecgram=1;

showExample=1;
useTrials=1:size(LFPbySweep,1);
UPthresh=alphaThresh;
useTrials=useTrials-useTrials(1)+1;

powerRatio=[];

if showExample
    dontShowSpecgram=0;
    getIntegralsForAlpha(LFPbySweep,exampleTrial,LFP_Fs,lowerBand,upperBand,1,UPthresh,smoothBin);
end
UPstates_LFP=cell(length(useTrials),1);
dontShowSpecgram=1;
powerRatio=cell(length(useTrials),1);
lowPower=cell(length(useTrials),1);
highPower=cell(length(useTrials),1);
for i=1:length(useTrials)
    [UPstates_LFP{i},powerRatio{i},lowPower{i},highPower{i}]=getIntegralsForAlpha(LFPbySweep,useTrials(i),LFP_Fs,lowerBand,upperBand,0,UPthresh,smoothBin);
end
UPstates=UPstates_LFP;
end

function [returnUPs,power_ratio,lowPower,highPower]=getIntegralsForAlpha(LFPbySweep,trials,LFP_Fs,lowerBand,upperBand,showExample,UP_thresh,smoothBin)

useRatio=0; % else use low power

% Get the parameters used for making LFP_specgram
[gb_params,freq,gabor]=gabor_morlet_config(LFP_Fs,[],[]);

[p,LFP_specgram]=makeWaveletSpecgram(LFPbySweep,trials,LFP_Fs);

% For each time point, plot integral (area under power spectrum)
% of frequencies greater than divide over the integral of frequencies 
% less than divide
% Use UP_thresh as threshold to identify UP states from this power ratio
% analysis
band1_lower=lowerBand(1); % band1 is lowerBand
band1_upper=lowerBand(2);
band2_lower=upperBand(1); % band2 is upperBand
band2_upper=upperBand(2);
divide_ind1_lower=0;
divide_ind1_upper=0;
divide_ind2_lower=0;
divide_ind2_upper=0;
for i=1:length(freq)
    if freq(i)>=band1_lower
        divide_ind1_lower=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band1_upper
        divide_ind1_upper=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band2_lower
        divide_ind2_lower=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band2_upper
        divide_ind2_upper=i;
        break
    end
end

high_integral=zeros(size(LFP_specgram,2),1);
size_low_int=divide_ind1_upper-divide_ind1_lower+1; % band1 is lowerBand
size_high_int=divide_ind2_upper-divide_ind2_lower+1; % band2 is upperBand
low_integral=zeros(size(LFP_specgram,2),1);
for i=1:size(LFP_specgram,2)
    low_area=sum(LFP_specgram(divide_ind1_lower:divide_ind1_upper,i)); % band1 is lowerBand
    high_area=sum(LFP_specgram(divide_ind2_lower:divide_ind2_upper,i)); % band2 is upperBand
    high_integral(i)=high_area;
    low_integral(i)=low_area;
end

% Normalize high and low integral by number of frequencies included
% power_ratio=(high_integral/size_high_int)./(low_integral/size_low_int);
% power_ratio=high_integral./low_integral;
power_ratio=low_integral./high_integral;
lowPower=low_integral;
highPower=high_integral;
if useRatio==0
    power_ratio=lowPower;
end
power_ratio=smooth(power_ratio,smoothBin);

times=0:(1/LFP_Fs):(1/LFP_Fs)*(size(LFPbySweep,2)-1);
if showExample
    figure(); 
    subplot(3,1,1);
    plot(times,LFPbySweep(trials,:));
    subplot(3,1,2);
    [Nfreq Ntime] = size(LFP_specgram);
    TimeLo=0;
    TimeHi=5;
    x = linspace(TimeLo,TimeHi,Ntime);
    TimeStep=0.5;
    Xtick = TimeLo : TimeStep : TimeHi;
    FreqLo=1;
    FreqHi=100;
    logy = linspace(log(FreqLo),log(FreqHi),Nfreq);
    Y_ticks=[1 4 8 10 20 30 40 50 60 80 100];
    Ytick = Y_ticks;
    imagesc(x,logy,LFP_specgram);
    set(gca, ...
        'XTick', Xtick, ...
        'XTickLabel', arrayfun(@num2str, Xtick, 'UniformOutput', false), ...
        'Ydir','normal', ...
        'YTick', log(Ytick), ...
        'YTickLabel', arrayfun(@num2str, Ytick, 'UniformOutput', false));
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    subplot(3,1,3);
    plot(times',power_ratio);
    ti=sprintf('Power Ratio');
    title(ti);
    xlim([times(1) times(end)]);
    xlabel('Time (s)');
    ylabel('Power Ratio');
    % Draw UP threshold
    line([times(1) times(end)],[UP_thresh UP_thresh],'Color','r');
end

% Return UPs
in_UP=0;
UP_starts=[];
UP_ends=[];
for i=1:size(LFP_specgram,2)
    if in_UP==0
        if power_ratio(i)>=UP_thresh
            in_UP=1;
            UP_starts=[UP_starts; times(i)];
        end
    else
        if power_ratio(i)<UP_thresh
            in_UP=0;
            UP_ends=[UP_ends; times(i-1)];
        end
    end
end
if in_UP==1
    UP_ends=[UP_ends; times(end)];
end
returnUPs=[UP_starts UP_ends];
end

function alphaStates=returnAlpha(xpoints,powerRatio,thresh)

alphaStates=cell(length(powerRatio),1);
for i=1:length(powerRatio)
    alphaStates{i}=returnStates(xpoints,powerRatio{i},thresh);
end
end

function UPs=returnStates(xpoints,ypoints,UP_thresh)
in_UP=0;
UP_starts=[];
UP_ends=[];
for i=1:length(xpoints)
    if in_UP==0
        if ypoints(i)>=UP_thresh
            in_UP=1;
            UP_starts=[UP_starts; xpoints(i)];
        end
    else
        if ypoints(i)<UP_thresh
            in_UP=0;
            UP_ends=[UP_ends; xpoints(i-1)];
        end
    end
end
if in_UP==1
    UP_ends=[UP_ends; xpoints(end)];
end
UPs=[UP_starts UP_ends];
end