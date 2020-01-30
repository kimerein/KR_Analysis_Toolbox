function dLGN_vs_ntsr1_analysis_for_psths(psth,V1psth,noThetaTrials,outputDir,V1dir,onlyUseTrials,addToName,stimWindow,sortByLGN)

% saveIndividualUnits=true;

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
% doF1analysis_saveAllTrialSpecs([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
% doF1analysis_saveAllTrialSpecs([],[],V1dir,[],V1psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);

a=load([outputDir '\noTheta_trialAv_temp_noLED.mat']);
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.noTheta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,1,stimWindow,sortByLGN);
firststep=find(noThetaTrials'==1 & ismember(psth.unitLED{1},ledOffVal));

% a=load([outputDir '\noTheta_noLED.mat']);
% temp=addAll(a.noTheta.allS);

currTrials=ismember(firststep,lowF1Trials);
dLGNpsth_noTheta_low_F1=cutPSTHtoOnlyTrials(psth,currTrials);
save([outputDir '\dLGNpsth_noTheta_low_F1' addToName '.mat'],'dLGNpsth_noTheta_low_F1');
save([outputDir '\trials_noTheta_low_F1' addToName '.mat'],'currTrials');

currTrials=ismember(firststep,highF1Trials);
dLGNpsth_noTheta_high_F1=cutPSTHtoOnlyTrials(psth,currTrials);
save([outputDir '\dLGNpsth_noTheta_high_F1' addToName '.mat'],'dLGNpsth_noTheta_high_F1');
save([outputDir '\trials_noTheta_high_F1' addToName '.mat'],'currTrials');


% V1
currTrials=ismember(firststep,lowF1Trials);
V1psth_noTheta_low_F1=cutPSTHtoOnlyTrials(V1psth,currTrials);
save([V1dir '\V1psth_noTheta_low_F1' addToName '.mat'],'V1psth_noTheta_low_F1');

currTrials=ismember(firststep,highF1Trials);
V1psth_noTheta_high_F1=cutPSTHtoOnlyTrials(V1psth,currTrials);
save([V1dir '\V1psth_noTheta_high_F1' addToName '.mat'],'V1psth_noTheta_high_F1');







a=load([outputDir '\theta_trialAv_temp_noLED.mat']);
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.theta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,0,stimWindow,sortByLGN);
firststep=find(noThetaTrials'==0 & ismember(psth.unitLED{1},ledOffVal));

% a=load([outputDir '\theta_noLED.mat']);
% temp=addAll(a.theta.allS);

currTrials=ismember(firststep,lowF1Trials);
dLGNpsth_theta_low_F1=cutPSTHtoOnlyTrials(psth,currTrials);
save([outputDir '\dLGNpsth_theta_low_F1' addToName '.mat'],'dLGNpsth_theta_low_F1');
save([outputDir '\trials_theta_low_F1' addToName '.mat'],'currTrials');

currTrials=ismember(firststep,highF1Trials);
dLGNpsth_theta_high_F1=cutPSTHtoOnlyTrials(psth,currTrials);
save([outputDir '\dLGNpsth_theta_high_F1' addToName '.mat'],'dLGNpsth_theta_high_F1');
save([outputDir '\trials_theta_high_F1' addToName '.mat'],'currTrials');


% V1
currTrials=ismember(firststep,lowF1Trials);
V1psth_theta_low_F1=cutPSTHtoOnlyTrials(V1psth,currTrials);
save([V1dir '\V1psth_theta_low_F1' addToName '.mat'],'V1psth_theta_low_F1');

currTrials=ismember(firststep,highF1Trials);
V1psth_theta_high_F1=cutPSTHtoOnlyTrials(V1psth,currTrials);
save([V1dir '\V1psth_theta_high_F1' addToName '.mat'],'V1psth_theta_high_F1');

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
