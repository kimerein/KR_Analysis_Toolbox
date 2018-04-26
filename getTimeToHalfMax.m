function halfTimes=getTimeToHalfMax(subData,xstep,fracRemaining,delay)

% subData from normalizeFR_turningOff_acrossChs
thresh=fracRemaining*200;
halfTimes=zeros(16,1);
for i=1:16
    a=find(subData(i,floor(delay/xstep):end)<thresh,1,'first');
    if isempty(a)
        halfTimes(i)=length(subData(i,floor(delay/xstep):end));
    else
        halfTimes(i)=a;
    end
end
halfTimes=halfTimes.*xstep;