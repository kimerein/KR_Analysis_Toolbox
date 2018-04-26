function result=pulseTriggeredLFP(pulseOnsetTime,LFPbySweep,LFP_Fs,includeTrials,color,fractionOfTotalTrial)
% Currently assumes 1 pulse per trial

pulseTrigLFPs=zeros(length(includeTrials),size(LFPbySweep,2));
midPointInd=floor(size(LFPbySweep,2)/2);
halfWindow=midPointInd-1;

alignInd=floor(pulseOnsetTime*LFP_Fs)+1;
upperInd=alignInd+halfWindow;
lowerInd=alignInd-halfWindow;
if upperInd>size(LFPbySweep,2)
    upperInd=size(LFPbySweep,2);
end
if lowerInd<1
    lowerInd=1;
end

result=mean(LFPbySweep(