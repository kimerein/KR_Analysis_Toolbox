function wrapperForF1analysis(psth,trialDuration,stimBySweep,outputDir,noThetaTrials)

stimConds=psth.unitStimcond{1};
stimVal=unique(stimConds);

% Align psth to vis stim
times=linspace(0,trialDuration,size(stimBySweep,2));
psth=alignPSTHtoStim(psth,stimBySweep,stimConds,stimVal,times,[]);

% Do F1 analysis
a=unique(psth.unitTrials{1});
for i=1:length(stimVal)
    uses{i}={[stimVal(i)]};
    uses_tri{i}={a};
end
disp('uses');
disp(uses);
disp('uses_tri');
disp(uses_tri);

led=unique(psth.unitLED{1});
ledOffVal=0;
ledOnVal=led(led>=0.01);

doF1analysis([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
