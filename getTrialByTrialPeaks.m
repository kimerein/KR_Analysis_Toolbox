function [ledOff_p,ledOn_p,peaks_ledOff_mean,peaks_ledOn_mean]=getTrialByTrialPeaks(psth,useUnits,binsizeForMax,noThetaTrials,isNoTheta)

if ~isempty(useUnits)
    psth.psths=psth.psths(useUnits);
end

default_bin=0.01;
% ntimes=floor(binsize/default_bin);
formax_ntimes=floor(binsizeForMax/default_bin);

ledOff=0;
ledOn=5.05;

l=psth.unitLED{1};
t=psth.t;
edgesForPeaks=t(1):binsizeForMax:t(end);

ledOff_p=cell(1,length(psth.psths));
ledOn_p=cell(1,length(psth.psths));
peaks_ledOff_mean=nan(length(psth.psths),length(edgesForPeaks)-1);
peaks_ledOn_mean=nan(length(psth.psths),length(edgesForPeaks)-1);

for i=1:length(psth.psths)
    p=psth.psths{i};
    ledOn_peaks=[];
    ledOff_peaks=[];
    for j=1:size(p,1)
        if noThetaTrials(j)~=isNoTheta
            continue
        else
            currTrial_peaks=nan(1,length(edgesForPeaks)-1);
            currTrial_psth=p(j,:);
            for k=1:length(edgesForPeaks)-1
                currTrial_peaks(k)=max(currTrial_psth(t>=edgesForPeaks(k) & t<edgesForPeaks(k+1)));
            end
        end
        if l(j)==0
            ledOff_peaks=[ledOff_peaks; currTrial_peaks];
        else
            ledOn_peaks=[ledOn_peaks; currTrial_peaks];
        end
    end
    ledOff_p{i}=ledOff_peaks;
    ledOn_p{i}=ledOn_peaks;
    peaks_ledOff_mean(i,:)=nanmean(ledOff_peaks,1);
    peaks_ledOn_mean(i,:)=nanmean(ledOn_peaks,1);
end


figure();
plot(linspace(0,14.5,size(peaks_ledOff_mean,2)),nanmean(peaks_ledOff_mean,1),'Color','k');
hold on;
plot(linspace(0,14.5,size(peaks_ledOn_mean,2)),nanmean(peaks_ledOn_mean,1),'Color','r');

            
