function [spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,startTime,stopTime,includeTrials,color,fractionOfTotalTrial)
% Calculates the spike-triggered LFP
% 
% spikes is a spike structure containing the spikes of interest
% 
% LFPbySweep is an array of LFP data
% Each row is a sweep; each column is a sample
% 
% LFP_Fs is the sampling rate of the LFP
% 
% only spikes with spike times between startTime and stopTime will be
% included in the spike-triggered average

% Get LFP around each spike
% Due to spikes near beginning and end of a sweep, will be some LFP blocks
% that are longer or shorter than the average
% Pad the shorter LFP blocks with NaN to allow averaging of the LFP
% Use a special average (specialMean) that ignores NaN
%
% fractionOfTotalTrial is the fraction of total sweep time to display
% surrounding spike
% e.g., if fractionOfTotalTrial = 1, will show whole trial

% Get current figure and current axes
hax = gca;

spikeTrigLFPs=zeros(length(spikes.spiketimes),size(LFPbySweep,2));
midPointInd=floor(size(LFPbySweep,2)/2);
halfWindow=midPointInd-1;
spikeN=1;
for i=1:max(spikes.trials)
    if any(includeTrials==i)
        tempspikes=filtspikes(spikes,0,'trials',i);
        for j=1:length(tempspikes.spiketimes)
            if tempspikes.spiketimes(j)>=startTime && tempspikes.spiketimes(j)<=stopTime
                LFPind=floor(tempspikes.spiketimes(j)/(1/LFP_Fs))+1;
                upperInd=LFPind+halfWindow;
                lowerInd=LFPind-halfWindow;
                if upperInd>size(LFPbySweep,2)
                    upperInd=size(LFPbySweep,2);
                end
                if lowerInd<1
                    lowerInd=1;
                end
                STLFP=LFPbySweep(i,lowerInd:upperInd);
                padLower=halfWindow-(LFPind-lowerInd);
                padUpper=halfWindow-(upperInd-LFPind)+1;
                paddedSTLFP=[NaN(1,padLower) STLFP NaN(1,padUpper)];
                spikeTrigLFPs(spikeN,:)=paddedSTLFP;
                spikeN=spikeN+1;
            end
        end
    end
end
temp_spikeTrigLFPs=spikeTrigLFPs(1:spikeN-1,:);
clear spikeTrigLFPs
spikeTrigLFPs=temp_spikeTrigLFPs;

% Plot the spike-triggered LFP and return average
% spike-triggered LFP
spikeTrigOut=spikeTrigLFPs(:,midPointInd-floor(halfWindow/(1/fractionOfTotalTrial)):midPointInd+floor(halfWindow/(1/fractionOfTotalTrial)));
spikeTrigAv=specialMean(spikeTrigOut,1);
h=plot(((size(spikeTrigOut,2)*(1/LFP_Fs))/2)-fliplr((1:size(spikeTrigOut,2))*(1/LFP_Fs)),spikeTrigAv,'Color',color);
%title('Spike-Triggered LFP');
% ylabel('LFP');
% xlabel('Time from Spike (s)');

function data_mean=specialMean(data, dim)
if dim==1
    data_sum=zeros(1,size(data,2));
    data_count=ones(1,size(data,2));
    for i=1:size(data,2)
        dataseg=data(:,i);
        if isempty((dataseg(~isnan(dataseg))))
            data_sum(i)=0;
            data_count(i)=1;
        else
            data_sum(i)=sum(dataseg(~isnan(dataseg)));
            data_count(i)=length(dataseg(~isnan(dataseg)));
        end
    end
    data_mean=data_sum./data_count;
elseif dim==2
    % need to fill in, but works without it
end

    