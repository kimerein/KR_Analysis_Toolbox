function [xpoints,ypoints]=psthWithOnlySpiketimes(spiketimes,windowSize,duration,Fs)

timeSteps=1/Fs;

for i=1:length(0:windowSize:duration)
    if i*windowSize>duration
        spikesInThisBin=sum((spiketimes>=(i-1)*windowSize) & (spiketimes<duration));
        xpoints(i)=(i-1)*windowSize;
        ypoints(i)=spikesInThisBin/(duration-((i-1)*windowSize));
    else
        spikesInThisBin=sum((spiketimes>=(i-1)*windowSize) & (spiketimes<i*windowSize));
        xpoints(i)=(i-1)*windowSize;
        ypoints(i)=spikesInThisBin/((i*windowSize)-((i-1)*windowSize));
    end
end

figure();
plot(xpoints,ypoints);