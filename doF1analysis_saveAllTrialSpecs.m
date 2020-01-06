function doF1analysis_saveAllTrialSpecs(expt,useFileInd,outputDir,LFP_Fs,dLGNpsth,usel_noLED,usel_LED,LFPdata,uses,uses_tri,noThetaTrials)

%uses=unique(dLGNpsth.unitStimcond{1});
usealls=0:10000;
% usealls=[1 1.05 2 2.05 4 4.05 6 6.05 8 8.05 10 10.05 12 12.05 14 14.05 16 16.05 18 18.05 20 20.05 30 30.05 40 40.05 50 50.05 60 60.05];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([outputDir '\' 'noThetaTrials.mat'],'noThetaTrials');
save([outputDir '\' 'dLGNpsth.mat'],'dLGNpsth');

disp('3');
tri=dLGNpsth.unitTrials{1};
if length(noThetaTrials)~=length(tri)
    disp('Length of noThetaTrials does not match number of trials in dLGNpsth');
    return
end
filt_psth=filtPSTH(dLGNpsth,tri(noThetaTrials'));
filt_psth_theta=filtPSTH(dLGNpsth,tri(noThetaTrials'~=1));

% No LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_noLED;

for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(single(triforsub),single(subcurrtri)) & ismember(single(sforsub),single(subcurrs))));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(single(triforsub),single(subcurrtri)) & ismember(single(sforsub),single(subcurrs))));
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
    end    
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel);
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_noLED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_noLED.mat'],'theta_trialAv_temp');
noTheta_trialAv=avResponses(noTheta_trialAv_temp);
theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_noLED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
save([outputDir '\' 'noTheta_noLED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
save([outputDir '\' 'theta_noLED.mat'],'theta');

disp('4');
% With LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usel=usel_LED;
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    allcurr_psth=[];
    allcurr_psth_theta=[];
    for j=1:length(currs)
        subcurrs=currs{j};
        subcurrtri=currtri{j};
        
        sforsub=filt_psth.unitStimcond{1};
        triforsub=filt_psth.unitTrials{1};
        curr_psth=filtPSTH(filt_psth,triforsub(ismember(single(triforsub),single(subcurrtri)) & ismember(single(sforsub),single(subcurrs))));
        
        sforsub=filt_psth_theta.unitStimcond{1};
        triforsub=filt_psth_theta.unitTrials{1};
        curr_psth_theta=filtPSTH(filt_psth_theta,triforsub(ismember(single(triforsub),single(subcurrtri)) & ismember(single(sforsub),single(subcurrs))));
        
        allcurr_psth=concatPSTHs(allcurr_psth,curr_psth);
        allcurr_psth_theta=concatPSTHs(allcurr_psth_theta,curr_psth_theta);
    end    
    [noTheta_trialAv_temp(i).allS,noTheta_trialAv_temp(i).HFa,noTheta_trialAv_temp(i).LFa,noTheta_trialAv_temp(i).F1amp,noTheta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth,ones(length(allcurr_psth.psths),1),zeros(length(allcurr_psth.psths),1),0,usealls,usel);
    [theta_trialAv_temp(i).allS,theta_trialAv_temp(i).HFa,theta_trialAv_temp(i).LFa,theta_trialAv_temp(i).F1amp,theta_trialAv_temp(i).allpower]=getHFandLFalphaResponses(allcurr_psth_theta,ones(length(allcurr_psth_theta.psths),1),zeros(length(allcurr_psth_theta.psths),1),0,usealls,usel);
end
save([outputDir '\' 'noTheta_trialAv_temp_LED.mat'],'noTheta_trialAv_temp');
save([outputDir '\' 'theta_trialAv_temp_LED.mat'],'theta_trialAv_temp');
noTheta_trialAv=avResponses(noTheta_trialAv_temp);
theta_trialAv=avResponses(theta_trialAv_temp);
save([outputDir '\' 'noTheta_trialAv_LED.mat'],'noTheta_trialAv');
save([outputDir '\' 'theta_trialAv_LED.mat'],'theta_trialAv');

[noTheta.allS,noTheta.HFa,noTheta.LFa,noTheta.F1amp,noTheta.allpower]=getHFandLFalphaResponses(filt_psth,ones(length(filt_psth.psths),1),zeros(length(filt_psth.psths),1),0,usealls,usel);
save([outputDir '\' 'noTheta_LED.mat'],'noTheta');
[theta.allS,theta.HFa,theta.LFa,theta.F1amp,theta.allpower]=getHFandLFalphaResponses(filt_psth_theta,ones(length(filt_psth_theta.psths),1),zeros(length(filt_psth_theta.psths),1),0,usealls,usel);
save([outputDir '\' 'theta_LED.mat'],'theta');

end

function outresponse=avResponses(response)

HFa=zeros(size(response(1).HFa));
LFa=zeros(size(response(1).LFa));
F1amp=zeros(size(response(1).F1amp));
allpower=zeros(size(response(1).allpower));
for i=1:length(response)
    if ~isnan(response(i).HFa)
        HFa=HFa+response(i).HFa;
    end
    if ~isnan(response(i).LFa)
        LFa=LFa+response(i).LFa;
    end
    if ~isnan(response(i).F1amp)
        F1amp=F1amp+response(i).F1amp;
    end
    if ~isnan(response(i).allpower)
        allpower=allpower+response(i).allpower;
    end
end
HFa=HFa./length(response);
LFa=LFa./length(response);
F1amp=F1amp./length(response);
allpower=allpower./length(response);
outresponse.HFa=HFa;
outresponse.LFa=LFa;
outresponse.F1amp=F1amp;
outresponse.allpower=allpower;

end


function newPSTH=concatPSTHs(struct1,struct2)

if isempty(struct1)
    newPSTH=struct2;
    return
elseif isempty(struct2)
    newPSTH=struct1;
    return
else
    if length(struct1.t)~=length(struct2.t)
        disp('t in psths does not match');
        return
    end
    newPSTH.t=struct1.t;
    newPSTH.psths=cell(length(struct1.psths),1);
    for i=1:length(struct1.psths)
        newPSTH.psths{i}=[struct1.psths{i}; struct2.psths{i}];
        newPSTH.unitTrials{i}=[struct1.unitTrials{i} struct2.unitTrials{i}];
        newPSTH.unitStimcond{i}=[struct1.unitStimcond{i} struct2.unitStimcond{i}];
        newPSTH.unitLED{i}=[struct1.unitLED{i} struct2.unitLED{i}];
    end
end

end


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
params.trialave=0;

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
        [S,t,f]=mtspecgramc(p,movingwin,params);
    else
        [S,t,f]=mtspecgramc(p(ismember(single(l),single(usel)) & ismember(single(s),single(uses)),:)',movingwin,params);
    end
    allS.S{i}=S;
    allS.t=t;
    allS.f=f;
    S=nanmean(S,3);
    S=reshape(S,size(S,1),size(S,2));
    HFa(i,:)=nanmean(S(:,f>=HFrange(1) & f<=HFrange(2)),2)';
    LFa(i,:)=nanmean(S(:,f>=LFrange(1) & f<=LFrange(2)),2)';
    F1amp(i,:)=nanmean(S(:,f>=F1range(1) & f<=F1range(2)),2)';
    allpower(i,:)=nanmean(S(:,f>=allrange(1) & f<=allrange(2)),2)';
end

end