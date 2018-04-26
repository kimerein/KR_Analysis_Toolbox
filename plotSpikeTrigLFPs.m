function f=plotSpikeTrigLFPs(spikes,LFPbySweep,bandPassedLFP,LFP_Fs,params,redTrials,blueTrials)

f=[];

f1=figure;
f=[f; f1];
% Spike-triggered LFP for all spikes
axesmatrix(6,3,1);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,0,params.totalTrialLength,1:max(spikes.trials),'k',0.1);
% Band-passed spike-triggered LFP for all spikes
axesmatrix(6,3,4);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,0,params.totalTrialLength,1:max(spikes.trials),'k',0.1);

% Spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,2);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,1:max(spikes.trials),'k',0.1);
% Band-passed spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,5);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,1:max(spikes.trials),'k',0.1);

% Spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,3);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,1:max(spikes.trials),'k',0.1);
% Band-passed spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,1:max(spikes.trials),'k',0.1);



% RED TRIALS
% Spike-triggered LFP for all spikes
axesmatrix(6,3,1+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,0,params.totalTrialLength,redTrials,'r',0.1);
% Band-passed spike-triggered LFP for all spikes
axesmatrix(6,3,4+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,0,params.totalTrialLength,redTrials,'r',0.1);

% Spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,2+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,redTrials,'r',0.1);
% Band-passed spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,5+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,redTrials,'r',0.1);

% Spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,3+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,redTrials,'r',0.1);
% Band-passed spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,6+6);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,redTrials,'r',0.1);



% BLUE TRIALS
% Spike-triggered LFP for all spikes
axesmatrix(6,3,1+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,0,params.totalTrialLength,blueTrials,'b',0.1);
% Band-passed spike-triggered LFP for all spikes
axesmatrix(6,3,4+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,0,params.totalTrialLength,blueTrials,'b',0.1);

% Spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,2+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,blueTrials,'b',0.1);
% Band-passed spike-triggered LFP for stimulus ON period spikes
axesmatrix(6,3,5+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart,params.ONstart+params.ONlength,blueTrials,'b',0.1);

% Spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,3+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,LFPbySweep,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,blueTrials,'b',0.1);
% Band-passed spike-triggered LFP for stimulus OFF period spikes
axesmatrix(6,3,6+12);
[spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(spikes,bandPassedLFP,LFP_Fs,params.ONstart+params.ONlength,params.totalTrialLength,blueTrials,'b',0.1);

nfigRows=length(params.Var1_name);
nfigCols=1;
% Plot spike-triggered LFPs for different orientations
if strcmp(params.Var1_name,'Orientation')
    f2=figure;
    f=[f; f2];
    % Spike-triggered LFP for orientation condition 1, etc.
    for i=1:length(params.Var1_values)
        someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+length(params.Var2_values)-1);
        axesmatrix(nfigRows,nfigCols,i);
        [spikeTrigOut,spikeTrigAv,h]=spikeTrigLFP(someSpikes,bandPassedLFP,LFP_Fs,0,params.totalTrialLength,someSpikes.sweeps.trials,'k',0.1);
        title(num2str(params.Var1_values(i)));
    end
end
    
        
        