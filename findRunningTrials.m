function [runningTrials]=findRunningTrials(velocity)

windowToAsses=[1.3 3];
diffForRunning=0.1;
trialDuration=5;

x=linspace(0,trialDuration,size(velocity,2));
runningTrials=zeros(size(velocity,1),1);
for i=1:size(velocity,1)
    currVec=smooth(velocity(i,:),300);
    subVec=currVec(x>=windowToAsses(1) & x<=windowToAsses(2));
    if max(subVec)-min(subVec)>diffForRunning
        runningTrials(i)=1;
    else
        runningTrials(i)=0;
    end
end
