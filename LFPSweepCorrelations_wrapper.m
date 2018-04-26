function LFPSweepCorrelations_wrapper(LFPbySweep,stimsForSweeps,LFP_Fs)

baselineStart=0;
baselineEnd=0.1;
stimStart=0.1;
stimEnd=1.1;
postStart=1.1;
postEnd=2;

stimsForSweeps(floor((stimsForSweeps-1)/6)+1==1)=1;
stimsForSweeps(floor((stimsForSweeps-1)/6)+1==2)=2;
stimsForSweeps(floor((stimsForSweeps-1)/6)+1==3)=3;
stimsForSweeps(floor((stimsForSweeps-1)/6)+1==4)=4;
stimsForSweeps(floor((stimsForSweeps-1)/6)+1==5)=5;

% timeVector=0:baselineEnd/(size(LFPbySweep(:,floor(baselineStart*LFP_Fs)+1:floor(baselineEnd*LFP_Fs)),2)-1):baselineEnd;
% LFPsweepCorrelations_withDistance(LFPbySweep(:,floor(baselineStart*LFP_Fs)+1:floor(baselineEnd*LFP_Fs)),stimsForSweeps,timeVector);

% timeVector=0:(stimEnd-stimStart)/(size(LFPbySweep(:,floor(stimStart*LFP_Fs)+1:floor(stimEnd*LFP_Fs)),2)-1):(stimEnd-stimStart);
% LFPsweepCorrelations_withDistance(LFPbySweep(:,floor(stimStart*LFP_Fs)+1:floor(stimEnd*LFP_Fs)),stimsForSweeps,timeVector);

timeVector=0:(postEnd-postStart)/(size(LFPbySweep(:,floor(postStart*LFP_Fs)+1:floor(postEnd*LFP_Fs)),2)-1):(postEnd-postStart);
LFPsweepCorrelations_withDistance(LFPbySweep(:,floor(postStart*LFP_Fs)+1:floor(postEnd*LFP_Fs)),stimsForSweeps,timeVector);