function [halfWidths,thinBools]=getSpikeWidths(spikes,Fs)

thresh=3.2*10^-4;

thinBools=zeros(1,size(spikes.waveforms,1));
halfWidths=zeros(1,size(spikes.waveforms,1));

for i=1:size(spikes.waveforms,1)
    [h,b]=classifyWaveformWidth(spikes,thresh,i,32000);
    thinBools(i)=b;
    halfWidths(i)=h;
    if mod(i,10000)==1
        disp(i);
    end
end

[n,xout]=hist(halfWidths,50);
figure();
bar(xout,n,1);