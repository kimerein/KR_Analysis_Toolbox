function [coh]=calcCoherenceContinuous(sineFreq,otherInput,stimWindow)

% stimWindow=[1.3 3];
Fs=1000;
useStimcond=1:128;
% useLED=sineFreq+0.050;

% params.tapers=[2.94 stimWindow(2)-stimWindow(1) 0];
% params.tapers=[0.1 30 0];
% params.tapers=[5 9];
params.tapers=[3 5];
params.Fs=Fs;
params.fpass=[0 70];
params.pad=0;
% params.err=[2 0.05];
% params.trialave=1;

stim.x=linspace(1,3,2/(1/Fs));
stim.y=sin(2*pi*floor(sineFreq).*stim.x);
stim.x=stim.x(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
stim.y=stim.y(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));

data1=stim.y';

% [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(data1,data2',params);
[C,phi,S12,S1,S2,f]=coherencycpb(data1,otherInput',params);
coh=mean(C(f>=sineFreq-0.4 & f<=sineFreq+0.4));
 
 