% function [allS,HFa,LFa,F1amp,allpower]=getHFandLFalphaResponses(psth1,likes1,sigCells,avBeforeSpec,uses,usel,currF1)
function [allS,HFa,LFa,F1amp,allpower]=getHFandLFalphaResponses(psth1,likes1,sigCells,avBeforeSpec,uses,usel)

allpower=[];
allS=[];
HFa=[];
LFa=[];
F1amp=[];

% usel=[0];
% uses=[1];
% HFrange=[11.5 12.5];
HFrange=[11.5 20];
LFrange=[4 6];
F1range=[2.5 3.5];
% F1range=[currF1-0.5 currF1+0.5];
% F1range=[currF1-0.75 currF1+0.75];
% if currF1==1
%     F1range=[1 2];
% end
% F1range=[1 10];
allrange=[10 99.5];
% avBeforeSpec=1;

movingwin=[1 0.05]; % v1
% movingwin=[0.165 0.025];
% movingwin=2*[0.165 0.025];
params.tapers=[5 6]; % v1
% params.tapers=[0.9 2 0];
% params.tapers=[3 5];
params.Fs=1/(psth1.t(2)-psth1.t(1));
% params.fpass=[1 30];
params.fpass=[1 100];
params.trialave=1;

% useUnits=find(likes1==1 & sigCells<=0.05);
useUnits=find(likes1==1 & sigCells<=1);
allS.t=[];
allS.f=[];
allS.S=cell(1,length(useUnits));
for i=1:length(useUnits)
    p=psth1.psths{useUnits(i)};
    l=psth1.unitLED{useUnits(i)};
    s=psth1.unitStimcond{useUnits(i)};
    if avBeforeSpec==1 
        p=nanmean(p(ismember(single(l),single(usel)) & ismember(single(s),single(uses)),:),1)';
        [S,t,f]=mtspecgrampb(p,movingwin,params);
    else
        [S,t,f]=mtspecgrampb(p(ismember(single(l),single(usel)) & ismember(single(s),single(uses)),:)',movingwin,params);
    end
    allS.S{i}=S;
    allS.t=t;
    allS.f=f;
    HFa(i,:)=nanmean(S(:,f>=HFrange(1) & f<=HFrange(2)),2)';
    LFa(i,:)=nanmean(S(:,f>=LFrange(1) & f<=LFrange(2)),2)';
    F1amp(i,:)=nanmean(S(:,f>=F1range(1) & f<=F1range(2)),2)';
    allpower(i,:)=nanmean(S(:,f>=allrange(1) & f<=allrange(2)),2)';
end

