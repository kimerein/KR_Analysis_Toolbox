function noThetaTrials=getBrainStateFromNoLED(thetaDiff,thetaThresh)

useWindow1=[0 3];
useWindow2=[9.5 14.5];
trialDuration=14.5;

times=linspace(0,trialDuration,size(thetaDiff,2));
trialBytrial=nanmean(thetaDiff(:,(times>=useWindow1(1) & times<=useWindow1(2)) | (times>=useWindow2(1) & times<=useWindow2(2))),2);

% Was THIS but note that I accidentally switched!
noThetaTrials=trialBytrial>=thetaThresh;

% noThetaTrials=trialBytrial<thetaThresh;