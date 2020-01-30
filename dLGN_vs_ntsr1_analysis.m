function dLGN_vs_ntsr1_analysis(psth,V1psth,noThetaTrials,outputDir,V1dir,onlyUseTrials,addToName,stimWindow,sortByLGN)

saveIndividualUnits=true;
saveFieldByField=true; % if files too big to save in one .mat

% outputDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Simultaneous dLGN and Ntsr1\Mawake401\dLGN\specs';
% V1dir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\FF_manuscript\Simultaneous dLGN and Ntsr1\Mawake401\V1 Ntsr1\specs';

s=unique(psth.unitStimcond{1});
clear uses; clear uses_tri;
a=psth.unitTrials{1};
if ~isempty(onlyUseTrials)
    a=a(ismember(a,onlyUseTrials));
    noThetaTrials=noThetaTrials(onlyUseTrials);
    psth=cutPSTHtoOnlyTrials(psth,onlyUseTrials);
    V1psth=cutPSTHtoOnlyTrials(V1psth,onlyUseTrials);
end
for i=1:length(s)
    uses{i}={[s(i)]};
    uses_tri{i}={a};
end
% ledOffVal=[0 0.05];
ledOffVal=[0];
led=unique(psth.unitLED{1});
ledOnVal=led(~ismember(led,ledOffVal));
ledOnVal=100; % to save time

% Comment out if already done
% doF1analysis_saveAllTrialSpecs([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials,saveFieldByField);
% doF1analysis_saveAllTrialSpecs([],[],V1dir,[],V1psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials,saveFieldByField);

a=loadVar(outputDir,'noTheta_trialAv_temp_noLED',saveFieldByField,'noTheta_trialAv_temp');
t=a.noTheta_trialAv_temp(1).allS.t;
f=a.noTheta_trialAv_temp(1).allS.f;
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.noTheta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,1,stimWindow,sortByLGN);
firststep=find(noThetaTrials'==1 & ismember(psth.unitLED{1},ledOffVal));

a=loadVar(outputDir,'noTheta_noLED',saveFieldByField,'noTheta');
if saveFieldByField==true
    a.noTheta.allS.S=loadVar(outputDir,'noTheta_noLED_allS_S',saveFieldByField,'S');
    a.noTheta.allS.t=t;
    a.noTheta.allS.f=f;
end
temp=addAll(a.noTheta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeLow'); title('No theta low F1');
save([outputDir '\noTheta_lowF1_dLGN' addToName '.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeHigh'); title('No theta high F1');
save([outputDir '\noTheta_highF1_dLGN' addToName '.mat'],'takeHigh');

% V1
a=loadVar(V1dir,'noTheta_noLED',saveFieldByField,'noTheta');
if saveFieldByField==true
    a.noTheta.allS.S=loadVar(V1dir,'noTheta_noLED_allS_S',saveFieldByField,'temp');
    b=loadVar(V1dir,'noTheta_trialAv_temp_noLED',saveFieldByField,'noTheta_trialAv_temp');
    a.noTheta.allS.t=b.noTheta_trialAv_temp(1).allS.t;
    a.noTheta.allS.f=b.noTheta_trialAv_temp(1).allS.f;
end
temp=addAll(a.noTheta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeLow'); title('V1 units -- No theta low F1');
save([V1dir '\noTheta_lowF1_Ntsr1' addToName '.mat'],'takeLow');
if saveIndividualUnits==true
    takeLow=[];
    for i=1:length(a.noTheta.allS.S)
        tempInd=a.noTheta.allS.S{i};
        takeLow=cat(3,takeLow,reshape(nanmean(tempInd(:,:,ismember(firststep,lowF1Trials)),3),size(tempInd,1),size(tempInd,2)));
    end
    save([V1dir '\noTheta_lowF1_Ntsr1_individUnits' addToName '.mat'],'takeLow');
end

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeHigh'); title('V1 units -- No theta high F1');
save([V1dir '\noTheta_highF1_Ntsr1' addToName '.mat'],'takeHigh');
if saveIndividualUnits==true
    takeHigh=[];
    for i=1:length(a.noTheta.allS.S)
        tempInd=a.noTheta.allS.S{i};
        takeHigh=cat(3,takeHigh,reshape(nanmean(tempInd(:,:,ismember(firststep,highF1Trials)),3),size(tempInd,1),size(tempInd,2)));
    end
    save([V1dir '\noTheta_highF1_Ntsr1_individUnits' addToName '.mat'],'takeHigh');
end




a=loadVar(outputDir,'theta_trialAv_temp_noLED',saveFieldByField,'theta_trialAv_temp');
t=a.theta_trialAv_temp(1).allS.t;
f=a.theta_trialAv_temp(1).allS.f;
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.theta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,0,stimWindow,sortByLGN);
firststep=find(noThetaTrials'==0 & ismember(psth.unitLED{1},ledOffVal));

a=loadVar(outputDir,'theta_noLED',saveFieldByField,'theta');
if saveFieldByField==true
    a.theta.allS.S=loadVar(outputDir,'theta_noLED_allS_S',saveFieldByField,'temp');
    a.theta.allS.t=t;
    a.theta.allS.f=f;
end
temp=addAll(a.theta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeLow'); title('Theta low F1');
save([outputDir '\theta_lowF1_dLGN' addToName '.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeHigh'); title('Theta high F1');
save([outputDir '\theta_highF1_dLGN' addToName '.mat'],'takeHigh');

% V1
a=loadVar(V1dir,'theta_noLED',saveFieldByField,'theta');
if saveFieldByField==true
    a.theta.allS.S=loadVar(V1dir,'theta_noLED_allS_S',saveFieldByField,'temp');
    b=loadVar(V1dir,'theta_trialAv_temp_noLED',saveFieldByField,'theta_trialAv_temp');
    a.theta.allS.t=b.theta_trialAv_temp(1).allS.t;
    a.theta.allS.f=b.theta_trialAv_temp(1).allS.f;
end
temp=addAll(a.theta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeLow'); title('V1 units -- Theta low F1');
save([V1dir '\theta_lowF1_Ntsr1' addToName '.mat'],'takeLow');
if saveIndividualUnits==true
    takeLow=[];
    for i=1:length(a.theta.allS.S)
        tempInd=a.theta.allS.S{i};
        takeLow=cat(3,takeLow,reshape(nanmean(tempInd(:,:,ismember(firststep,lowF1Trials)),3),size(tempInd,1),size(tempInd,2)));
    end
    save([V1dir '\theta_lowF1_Ntsr1_individUnits' addToName '.mat'],'takeLow');
end

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeHigh'); title('V1 units -- Theta high F1');
save([V1dir '\theta_highF1_Ntsr1' addToName '.mat'],'takeHigh');
if saveIndividualUnits==true
    takeHigh=[];
    for i=1:length(a.theta.allS.S)
        tempInd=a.theta.allS.S{i};
        takeHigh=cat(3,takeHigh,reshape(nanmean(tempInd(:,:,ismember(firststep,highF1Trials)),3),size(tempInd,1),size(tempInd,2)));
    end
    save([V1dir '\theta_highF1_Ntsr1_individUnits' addToName '.mat'],'takeHigh');
end

end

function a=loadVar(varargin)

if length(varargin)==3
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=[];
elseif length(varargin)==4
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=varargin{4};
else
    error('expected 3 or 4 arguments to loadVar');
end

if loadFieldByField==true
    if isempty(inFileName)
        a.(varName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(varName))
            a=a.(varName);
        end
    else
        a.(inFileName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(inFileName))
            a=a.(inFileName);
        end
    end
else
    a=load([varDir '\' varName '.mat']);
end

end

function currpsth=cutPSTHtoOnlyTrials(currpsth,takeTri)

for j=1:length(currpsth.psths)
    temp=currpsth.psths{j};
    temp=temp(takeTri,:);
    currpsth.psths{j}=temp;
    
    temp=currpsth.unitTrials{j};
    temp=temp(takeTri);
    currpsth.unitTrials{j}=temp;
    
    temp=currpsth.unitStimcond{j};
    temp=temp(takeTri);
    currpsth.unitStimcond{j}=temp;
    
    temp=currpsth.unitLED{j};
    temp=temp(takeTri);
    currpsth.unitLED{j}=temp;
end
    
end

function temp=addAll(allS)

temp=[];
for i=1:length(allS.S)
    tmp=cat(4,temp,allS.S{i});
    temp=nansum(tmp,4);
end
temp=temp./length(allS.S);

end

function [lowF1Trials,highF1Trials]=getHighVsLowF1(specs_temp,psth,uses,uses_tri,noThetaTrials,isNoTheta,stimWindow,sortByLGN)

if sortByLGN==true
    [lowF1Trials,highF1Trials]=findTrialsWithHighOrLowF1(specs_temp,psth,uses,uses_tri,noThetaTrials,isNoTheta,stimWindow);
else
    [lowF1Trials,highF1Trials]=findTrialsWithHighOrLowNtsr1(specs_temp,psth,uses,uses_tri,noThetaTrials,isNoTheta,stimWindow);
end

end
