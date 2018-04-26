function [dtSteps,fracSupp_spontAc,noSupp_spontAc,supp_spontAc]=calculatingSilencingForAllWindows(spikes,waveforms,infoDOTdetect,halfWidths,windowOfInterest)

% Defining parameters for analysis
detectThresh=3.2*10^-4; % threshold for separating thin and thick spiking
                        % waveforms based on half-width-at-half-max
spontSil_LEDconds=[7]; % LED conds for silencing of spontaneous activity
spontNo_LEDconds=[2 4]; % LED conds for no LED during spontaneous activity
% dtSteps=-0.04:0.003:0.06; % incremental changes to spike detection threshold
dtSteps=0.003;
trialDuration=6;
psthBinSize=0.01;
Fs=32000;
% windowOfInterest=[3.3 3.5];
thickSpikes=true; % Use only spikes with waveform half-width-at-half-max >= detectThresh
nChannels=16;
eventChToWvfrmMapping=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
% eventChToWvfrmMapping=1:16;

% Get waveform mins
minWvfrms=min(waveforms,[],2);

% Calculating silencing
supp_spontAc=zeros(length(dtSteps),1);
noSupp_spontAc=zeros(length(dtSteps),1);
fracSupp_spontAc=zeros(length(dtSteps),1);
for i=1:length(dtSteps)
    disp(dtSteps(i));
    detectCutoffs=zeros(size(min(minWvfrms(:,:,:),[],3)));
    for j=1:nChannels
        th=infoDOTdetect.thresh(j)-dtSteps(i);
        detectCutoffs(min(minWvfrms(:,:,eventChToWvfrmMapping(j)),[],3)<th & infoDOTdetect.event_channel==j)=1;
    end
    if thickSpikes % Comment out for no sorting of thick and thin spikes
        takeTheseWvfrms=halfWidths>=detectThresh;
    else
        takeTheseWvfrms=halfWidths<detectThresh;
    end
%     takeTheseWvfrms=ones(size(spikes.trials)); % Uncomment for no sorting of thick and thin spikes
    numT=length(unique(spikes.trials((takeTheseWvfrms)' & detectCutoffs & ismember(spikes.led,spontSil_LEDconds)')));
    subSpikes=spikes.spiketimes((takeTheseWvfrms)' & detectCutoffs & ismember(spikes.led,spontSil_LEDconds)');
    [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,psthBinSize,trialDuration,Fs);
%     if i~=1 close all; end
    supp_spontAc(i)=mean(ypoints(xpoints>windowOfInterest(1) & xpoints<windowOfInterest(2)))/(numT*(windowOfInterest(2)-windowOfInterest(1)));
    numT=length(unique(spikes.trials((takeTheseWvfrms)' & detectCutoffs & ismember(spikes.led,spontNo_LEDconds)')));
    subSpikes=spikes.spiketimes((takeTheseWvfrms)' & detectCutoffs & ismember(spikes.led,spontNo_LEDconds)');
    [xpoints,ypoints]=psthWithOnlySpiketimes(subSpikes,psthBinSize,trialDuration,Fs);
%     if i~=1 close all; end
    noSupp_spontAc(i)=mean(ypoints(xpoints>windowOfInterest(1) & xpoints<windowOfInterest(2)))/(numT*0.2);
    fracSupp_spontAc(i)=supp_spontAc(i)./noSupp_spontAc(i);
end

% Make figure of fractional suppression as function of spike detection
% threshold
% figure();
% plot(dtSteps,fracSupp_spontAc);
% figure();
% plot(dtSteps,noSupp_spontAc,'Color','k');
% hold on;
% plot(dtSteps,supp_spontAc,'Color','r');

%STATS=testt(ns_evokedNo_thick-m_evokedNo_Base_thick,ns_evoked4Sil_thick-m_evoked4Sil_Base_thick,0,0.05,2);