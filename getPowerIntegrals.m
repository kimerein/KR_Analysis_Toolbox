function [returnUPs,power_ratio]=getPowerIntegrals(LFPbySweep,trials,LFP_Fs,lowerBand,upperBand,showExample,UP_thresh,spikes,MUAthresh)

% Get the parameters used for making LFP_specgram
[gb_params,freq,gabor]=gabor_morlet_config(LFP_Fs,[],[]);

[p,LFP_specgram]=makeWaveletSpecgram(LFPbySweep,trials,LFP_Fs);

if ~isempty(spikes)
%     a=unique(spikes.trials);
    a=unique(spikes.sweeps.trials);
    if length(a)~=size(LFPbySweep,1)
        disp('spikes # of trials and size of LFPbySweep do not match!');
    end
    % The following line also in find_UPstates_with_LFP -- must be the same
    [n,c,e,xpoints,ypoints]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'trials',a(trials)),250,0);
%     [n,c,e,xpoints,ypoints]=psth_wStdev_valuesOnly(filtspikes(spikes,0,'trials',trials),250,0);
end

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
power_ratio=smooth(power_ratio,1500);

times=0:(1/LFP_Fs):(1/LFP_Fs)*(size(LFPbySweep,2)-1);
if showExample
    figure(); 
    plot(times',power_ratio);
    ti=sprintf('Power Ratio');
    title(ti);
    xlim([times(1) times(end)]);
    xlabel('Time (s)');
    ylabel('Power Ratio');
    % Draw UP threshold
    line([times(1) times(end)],[UP_thresh UP_thresh],'Color','r');
    if ~isempty(spikes)
        figure(); 
        plot(xpoints,ypoints);
        line([xpoints(1) xpoints(end)],[MUAthresh MUAthresh],'Color','r');
    end
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
%             UP_ends=[UP_ends; times(i)];
            UP_ends=[UP_ends; times(i-1)];
        end
    end
end
if in_UP==1
    UP_ends=[UP_ends; times(end)];
end
returnUPs=[UP_starts UP_ends];