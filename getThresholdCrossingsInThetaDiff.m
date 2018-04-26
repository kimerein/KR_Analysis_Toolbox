function output=getThresholdCrossingsInThetaDiff(thetaDiff,thetaThresh,trialDuration)

nCrossings=zeros(1,size(thetaDiff,1));
for i=1:size(thetaDiff,1)
    curr=thetaDiff(i,:);
    curr=downSampAv(curr,3);
    nCrossings(i)=sum(curr(1:end-1)<thetaThresh & curr(2:end)>thetaThresh);
end

% Frequency of switches
% in Hz
totalSeconds=size(thetaDiff,1)*trialDuration;
freqOfSwitch=sum(nCrossings)/totalSeconds;
output.freqOfSwitch=freqOfSwitch;

% Fraction of time in theta
output.fracTheta=sum(thetaDiff(1:end)>=thetaThresh)./length(thetaDiff(1:end));

