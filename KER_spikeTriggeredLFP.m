function spikeTrigOut=KER_spikeTriggeredLFP(spikes,LFPbySweep,LFP_Fs,startTime,stopTime)

trialLength=size(LFPbySweep,2)*(1/LFP_Fs);
spikeTrigLFPs=zeros(length(spikes.spiketimes),size(LFPbySweep,2));
midPointInd=floor(size(LFPbySweep,2)/2);
halfWindow=midPointInd-1;

spikeN=1;
for i=1:max(spikes.trials)
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
temp_spikeTrigLFPs=spikeTrigLFPs(1:spikeN-1,:);
clear spikeTrigLFPs
spikeTrigLFPs=temp_spikeTrigLFPs;

figure;
spikeTrigOut=spikeTrigLFPs(:,midPointInd-floor(halfWindow/100):midPointInd+floor(halfWindow/100));
%plot((trialLength/2)-fliplr(0:trialLength/(size(LFPbySweep,2)-1):trialLength),specialMean(spikeTrigLFPs,1));
plot(((size(spikeTrigOut,2)*(1/LFP_Fs))/2)-fliplr((1:size(spikeTrigOut,2))*(1/LFP_Fs)),specialMean(spikeTrigOut,1));
title('Spike-Triggered LFP');
ylabel('LFP');
xlabel('Time from Spike (s)');
end

function data_mean=specialMean(data, dim)
% if dim==1
%     otherdim=2;
% else
%     otherdim=1;
% end
% data_sum=zeros(1,size(data,otherdim));
% data_count=ones(1,size(data,otherdim));
% for i=1:size(data,dim)
%     for j=1:size(data,otherdim)
%         if
% end

% if dim==1
%     data_sum=zeros(1,size(data,2));
%     data_count=zeros(1,size(data,2));
%     for i=1:size(data,2)
%         innerSum=0;
%         innerCount=0;
%         for j=1:size(data,1)
%             if ~isequal(data(j,i),NaN)
%                 innerSum=innerSum+data(j,i);
%                 innerCount=innerCount+1;
%             end
%         end
%         if innerCount==0
%             data_sum(i)=0;
%             data_count(i)=1;
%         else
%             data_sum(i)=innerSum;
%             data_count(i)=innerCount;
%         end
%     end
% end
% data_mean=data_sum./data_count;

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
end
end

    