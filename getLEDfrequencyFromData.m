function LED=getLEDfrequencyFromData(LEDbySweep,Fs)

rmpath(genpath('chronux_2_11'));

peakThresh=1.5;
minTimeDiff=0.01;
stimWindow=[1 3];

L=LEDbySweep{1};

times=0:1/Fs:(size(L,2)-1)*(1/Fs);

peak_counts=nan(size(L,1),1);
for i=1:size(L,1)
    curr=L(i,:);
    peak_counts(i)=length(findpeaks(curr(times>stimWindow(1) & times<stimWindow(1)+1),Fs,'MinPeakHeight',peakThresh,'MinPeakDistance',minTimeDiff,'MinPeakProminence',0.5));
end

LED=peak_counts;