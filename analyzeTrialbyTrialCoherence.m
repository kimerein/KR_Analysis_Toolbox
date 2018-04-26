function [psths1,psths2,cohereAns]=analyzeTrialbyTrialCoherence(spikes1,spikes2,psths1,psths2,tau)

binsize=1; % in ms
averagePSTHs=1;
doConv=0;
trialDuration=10.5; % in s
% params.tapers=[3 5];
params.tapers=[0.2857 10.5 0];
params.Fs=1000/binsize;
params.fpass=[1 70];
params.pad=0;
params.err=[2 0.05];
params.trialave=1;

t=unique(spikes1.trials);
t2=unique(spikes2.trials);
if any(~ismember(t,t2))
    disp('trials do not match in two spikes structs');
    return
end

if isempty(psths1)
    [~,~,~,x,~,psths1]=psth_wStd_trialByTrial(spikes1,binsize,0,trialDuration,length(t),t);
    disp('Done with psths1');
end
if isempty(psths2)
    [~,~,~,x,~,psths2]=psth_wStd_trialByTrial(spikes2,binsize,0,trialDuration,length(t),t);
    disp('Done with psths2');
end

if averagePSTHs==1
    set1=1:2:size(psths1,1);
    set2=2:2:size(psths1,1);
    p1_set1=nanmean(psths1(set1,:),1);
    p2_set1=nanmean(psths2(set1,:),1);
    p1_set2=nanmean(psths1(set2,:),1);
    p2_set2=nanmean(psths2(set2,:),1);
%     if doConv==1
%         [out,result]=predictCxResponse_fromThalResponse(
    [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencypb([p1_set1; p1_set2]',[p2_set1; p2_set2]',params);
else
    [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencypb(psths1',psths2',params);
end
cohereAns.C=C;
cohereAns.phi=phi;
cohereAns.S12=S12;
cohereAns.S1=S1;
cohereAns.S2=S2;
cohereAns.f=f;
cohereAns.zerosp=zerosp;
cohereAns.confC=confC;
cohereAns.phistd=phistd;
cohereAns.Cerr=Cerr;

figure(); 
hax=axes();
hl=semilogx(f,nanmean(C,2),'Color','k');
hold on;
semilogx(f,Cerr(1,:),'Color','k');
semilogx(f,Cerr(2,:),'Color','k'); 