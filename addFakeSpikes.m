function [outspikes,outassignsinfo]=addFakeSpikes(spikes,newassignsinfo)

% cellmode='frequency double';
% cellmode='follow';
% cellmode='followx2';
cellmode='chirp';
% stimOn=[1 3];
stimOn=[0.5 9.5];

if strcmp(cellmode,'chirp')
    bin=1;
    stim_x=linspace(0,3,3*(1000/bin));
    stim_y_temp=chirp(stim_x,1,3,30,'logarithmic')+0.5;
    stim_y=[stim_y_temp fliplr(stim_y_temp) stim_y_temp];
    [~,chirpPeakInds]=findpeaks(stim_y);
    chirpfreq=1*((30/1)^(1/3)).^stim_x;
    chirpfreq=[chirpfreq fliplr(chirpfreq) chirpfreq];
    stim_x=linspace(0,9,9*(1000/bin));
    chirpSpiketimes=stimOn(1)+stim_x(chirpPeakInds);
end   

    leds=unique(spikes.led);
trials=spikes.sweeps.trials;
trialLEDs=spikes.sweeps.led;
trialStimconds=spikes.sweeps.stimcond;
trialFileinds=spikes.sweeps.fileInd;
trialTriggers=spikes.sweeps.trigger;
newass=max(unique(spikes.assigns))+1;
outspikes=spikes;
disp('totaltrials');
disp(length(trials));
for i=1:length(trials)
    if mod(i,500)==0
        disp(i);
    end
    currfreq=floor(trialLEDs(i));
    currStimcond=trialStimconds(i);
    currFileind=trialFileinds(i);
    currTrigger=trialTriggers(i);
    doSpiketimes=[];
    if strcmp(cellmode,'chirp')==1
        doSpiketimes=single(chirpSpiketimes);
    end
    if strcmp(cellmode,'frequency double')==1
        newFreq=2*currfreq;
        period=1/newFreq;
        for j=1:10000
            newSpiketime=j*period-(period/2);
            doSpiketimes=[doSpiketimes single(newSpiketime)];
            if newSpiketime>stimOn(2)
                break
            end
        end
        doSpiketimes=doSpiketimes+stimOn(1);      
    end
    if strcmp(cellmode,'follow')==1
        newFreq=currfreq;
        period=1/newFreq;
        for j=1:10000
            newSpiketime=j*period-(period/2);
            doSpiketimes=[doSpiketimes single(newSpiketime)];
            if newSpiketime>stimOn(2)
                break
            end
        end
        doSpiketimes=doSpiketimes+stimOn(1);      
    end
    if strcmp(cellmode,'followx2')==1
        newFreq=currfreq;
        period=1/newFreq;
        for j=1:10000
            newSpiketime=j*period-(period/2);
            doSpiketimes=[doSpiketimes single(newSpiketime) single(newSpiketime)];
            if newSpiketime>stimOn(2)
                break
            end
        end
        doSpiketimes=doSpiketimes+stimOn(1);      
    end
    outspikes.trials(size(outspikes.trials,2)+1:size(outspikes.trials,2)+1+length(doSpiketimes)-1)=single(trials(i)*ones(1,length(doSpiketimes)));
    outspikes.waveforms(size(outspikes.waveforms,1)+1:size(outspikes.waveforms,1)+1+length(doSpiketimes)-1,:,:)=single(zeros([length(doSpiketimes) size(outspikes.waveforms,2) size(outspikes.waveforms,3)]));
    outspikes.spiketimes(size(outspikes.spiketimes,2)+1:size(outspikes.spiketimes,2)+1+length(doSpiketimes)-1)=single(doSpiketimes);
    outspikes.unwrapped_times(size(outspikes.unwrapped_times,2)+1:size(outspikes.unwrapped_times,2)+1+length(doSpiketimes)-1)=single(0*ones(1,length(doSpiketimes)));
    outspikes.fileInd(size(outspikes.fileInd,2)+1:size(outspikes.fileInd,2)+1+length(doSpiketimes)-1)=single(currFileind*ones(1,length(doSpiketimes)));
    outspikes.trigger(size(outspikes.trigger,2)+1:size(outspikes.trigger,2)+1+length(doSpiketimes)-1)=single(currTrigger*ones(1,length(doSpiketimes)));
    outspikes.time(size(outspikes.time,2)+1:size(outspikes.time,2)+1+length(doSpiketimes)-1)=single(0*ones(1,length(doSpiketimes)));
    outspikes.stimcond(size(outspikes.stimcond,2)+1:size(outspikes.stimcond,2)+1+length(doSpiketimes)-1)=single(currStimcond*ones(1,length(doSpiketimes)));
    outspikes.led(size(outspikes.led,2)+1:size(outspikes.led,2)+1+length(doSpiketimes)-1)=single(trialLEDs(i)*ones(1,length(doSpiketimes)));
    outspikes.assigns(size(outspikes.assigns,2)+1:size(outspikes.assigns,2)+1+length(doSpiketimes)-1)=single(newass*ones(1,length(doSpiketimes)));
    if isfield(outspikes,'temp')
        outspikes.temp(size(outspikes.temp,2)+1:size(outspikes.temp,2)+1+length(doSpiketimes)-1)=single(0*ones(1,length(doSpiketimes)));
    end
end

if ~isempty(newassignsinfo)
    outassignsinfo=newassignsinfo;
    currSize=length(newassignsinfo.trode);
    outassignsinfo.trode(currSize+1)=outassignsinfo.trode(currSize);
    outassignsinfo.event_channel(currSize+1,:)=outassignsinfo.event_channel(currSize,:);
    outassignsinfo.depthTrode(currSize+1)=outassignsinfo.depthTrode(currSize);
    outassignsinfo.calibrated_evCh(currSize+1)=outassignsinfo.calibrated_evCh(currSize);
    outassignsinfo.waveforms(currSize+1,:)=outassignsinfo.waveforms(currSize,:);
    outassignsinfo.waveformWidths(currSize+1)=outassignsinfo.waveformWidths(currSize);
    outassignsinfo.original_assigns(currSize+1)=newass;
    outassignsinfo.new_assigns(currSize+1)=newass;
else
    outassignsinfo=[];
end
