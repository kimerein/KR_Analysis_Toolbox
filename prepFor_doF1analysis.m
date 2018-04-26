function prepFor_doF1analysis(spikes,thetaDiff,thetaThresh,taketri,saveDir,noLED,LED,uses,uses_tri)

% noThetaTrials=getBrainStateFromNoLED(thetaDiff,thetaThresh);
noThetaTrials=logical(zeros(1,length(taketri)));

dLGNspikes_subtrials=filtspikes(spikes,0,'trials',taketri);
dLGNspikes_subtrials.trials=dLGNspikes_subtrials.trials-min(dLGNspikes_subtrials.trials)+1;
dLGNspikes_subtrials.sweeps.trials=dLGNspikes_subtrials.sweeps.trials-min(dLGNspikes_subtrials.sweeps.trials)+1;
[psth]=measureAllUnitsResponseProperties(dLGNspikes_subtrials,unique(dLGNspikes_subtrials.assigns),[0 14.5]);

doF1analysis([],[],saveDir,[],psth,noLED,LED,[],uses,uses_tri,noThetaTrials);