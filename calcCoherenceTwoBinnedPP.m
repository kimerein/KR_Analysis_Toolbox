function [coh,ph,tdelay,C,f]=calcCoherenceTwoBinnedPP(sineFreq,firstInput,otherInput)

% if ismember(sineFreq,[1 2 4 6])
%     sineFreq=2*sineFreq;
% end

stimWindow=[1.3 3];
Fs=1000;
doCrossCorrFirst=0;
offsetBy3ms=0; % doesn't matter
% useLED=sineFreq+0.050;

% params.tapers=[2.94 stimWindow(2)-stimWindow(1) 0];
% params.tapers=[0.1 30 0];
% params.tapers=[5 9];
% params.tapers=[3 5];
params.tapers=[15 18];
% params.tapers=[18 21];
params.Fs=Fs;
params.fpass=[1 70];
params.pad=0;
params.trialave=1;
% params.err=[2 0.05];
% params.trialave=1;

% stim.x=linspace(1,3,2/(1/Fs));
% stim.y=sin(2*pi*floor(sineFreq).*stim.x);
% stim.x=stim.x(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
% stim.y=stim.y(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
% 
% data1=stim.y';

% [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(data1,data2',params);

if doCrossCorrFirst==1
    stim.x=linspace(1,3,2*Fs);
    stim.y=sin(2*pi*floor(sineFreq).*stim.x);
    stim.y=stim.y(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
    stim.x=stim.x(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));    
    [cc,cc_lags]=xcorr(stim.y-mean(stim.y),firstInput-mean(firstInput),[],'coeff');
    new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
    useFirstInput=cc(new_ccX>0);
    stim.x=linspace(1,3,2*Fs);
    stim.y=sin(2*pi*floor(sineFreq).*stim.x);
    stim.y=stim.y(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));
    stim.x=stim.x(stim.x>=stimWindow(1) & stim.x<=stimWindow(2));    
    [cc,cc_lags]=xcorr(stim.y-mean(stim.y),otherInput-mean(otherInput),[],'coeff');
    new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
    useOtherInput=cc(new_ccX>0);
    firstInput=useFirstInput;
    if offsetBy3ms==1
        otherInput=[useOtherInput(4:end) zeros(1,3)];
    else
        otherInput=useOtherInput;
    end
end

coh=[];
ph=[];
tdelay=[];
[C,phi,S12,S1,S2,f]=coherencypb(firstInput',otherInput',params);
% coh=mean(C(f>=sineFreq-0.4 & f<=sineFreq+0.4));
% ph=mean(phi(f>=sineFreq-0.4 & f<=sineFreq+0.4));
% tdelay=(ph/(2*pi))*(1/sineFreq);
 