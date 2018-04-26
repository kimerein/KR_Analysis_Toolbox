function [halfWidth,boolThin]=classifyUnitWaveformWidth(waveform,thresh,Fs)
showFigs=0;

shift1=double(-waveform)-min(double(-waveform));
[peak,peakInd]=findpeaks(shift1,'SORTSTR','descend','NPEAKS',1);
halfAmp=(peak-shift1(1))/2;
shift2=shift1-halfAmp;
point1=find(shift2(peakInd:-1:1)<0,1,'first');
if isempty(point1)
    point1=1;
end
point1=peakInd-point1+1;
point2=find(shift2(peakInd:end)<0,1,'first');
if isempty(point2)
    point2=length(shift2(peakInd:end));
end
point2=peakInd+point2-1;
halfWidth=point2-point1;
halfWidth=halfWidth/Fs;
if halfWidth<thresh
    boolThin=true;
else
    boolThin=false;
end

if showFigs==1
    figure(); 
    plot(waveform);
    hold on; 
    scatter(halfWidth*Fs,0);
end

