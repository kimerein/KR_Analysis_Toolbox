function [halfWidths,thinBools]=getUnitsSpikeWidths(spikes)
% WRONG
thresh=3.2*10^-4;
Fs=32000;

a=unique(spikes.assigns);

thinBools=zeros(1,length(a));
halfWidths=zeros(1,length(a));

for i=1:length(a)
    subSpikes=filtspikes(spikes,0,'assigns',a(i));
    maxEvCh=round(mean(subSpikes.info.detect.event_channel));
    wvfrm=mean(subSpikes.waveforms(:,:,maxEvCh),1);
    [h,b]=classifyUnitWaveformWidth(wvfrm,thresh,Fs);
    thinBools(i)=b;
    halfWidths(i)=h;
    disp(i);
end

[n,xout]=hist(halfWidths,10);
figure();
bar(xout,n,1);