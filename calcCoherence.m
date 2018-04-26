function [coh,cerr1,cerr2]=calcCoherence(sineFreq,spikes,useLED)

stimWindow=[1.3 3];
Fs=32000;
useStimcond=1:128;
% useLED=sineFreq+0.050;

% params.tapers=[2.94 stimWindow(2)-stimWindow(1) 0];
% params.tapers=[0.1 30 0];
% params.tapers=[5 9];
params.tapers=[5 9];
params.Fs=Fs;
params.fpass=[0 70];
params.pad=0;
params.err=[2 0.05];
params.trialave=1;
params.fscorr=1;
params.t=stimWindow(2)-stimWindow(1); % in s

stim.x=linspace(1,3,2/(1/Fs));
stim.y=sin(2*pi*floor(sineFreq).*stim.x);
stim.x=stim.x(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
stim.y=stim.y(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));

ttt=unique(spikes.trials(ismember(spikes.stimcond,useStimcond) & ismember(floor(spikes.led),floor(useLED))));
st=cell(length(ttt),1);
for i=1:length(ttt)
    currTrial=ttt(i);
    st{i}=spikes.spiketimes(spikes.trials==currTrial)'-stimWindow(1);
end
data1=repmat(stim.y',1,length(ttt));
data2=struct('times',st);
% data1=stim.y';
% data2=spikes.spiketimes(ismember(spikes.stimcond,useStimcond) & ismember(floor(spikes.led),floor(useLED)));

% [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(data1,data2',params);
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(data1,data2',params);
coh=mean(C(f>=sineFreq-0.4 & f<=sineFreq+0.4));
cerr1=mean(Cerr(1,f>=sineFreq-0.4 & f<=sineFreq+0.4));
cerr2=mean(Cerr(2,f>=sineFreq-0.4 & f<=sineFreq+0.4));
 
 