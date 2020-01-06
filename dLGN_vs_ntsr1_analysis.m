function dLGN_vs_ntsr1_analysis(psth,V1psth,noThetaTrials,outputDir,V1dir,onlyUseTrials,stimWindow)

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
ledOffVal=[0 0.05];
led=unique(psth.unitLED{1});
ledOnVal=led(~ismember(led,ledOffVal));
ledOnVal=100; % to save time

% Comment out if already done
doF1analysis_saveAllTrialSpecs([],[],outputDir,[],psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);
doF1analysis_saveAllTrialSpecs([],[],V1dir,[],V1psth,ledOffVal,ledOnVal,[],uses,uses_tri,noThetaTrials);

a=load([outputDir '\noTheta_trialAv_temp_noLED.mat']);
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.noTheta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,1,stimWindow);
firststep=find(noThetaTrials'==1 & ismember(psth.unitLED{1},ledOffVal));

a=load([outputDir '\noTheta_noLED.mat']);
temp=addAll(a.noTheta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeLow'); title('No theta low F1');
save([outputDir '\noTheta_lowF1_dLGN.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeHigh'); title('No theta high F1');
save([outputDir '\noTheta_highF1_dLGN.mat'],'takeHigh');

% V1
a=load([V1dir '\noTheta_noLED.mat']);
temp=addAll(a.noTheta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeLow'); title('V1 units -- No theta low F1');
save([V1dir '\noTheta_lowF1_Ntsr1.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.noTheta.allS.t,a.noTheta.allS.f,takeHigh'); title('V1 units -- No theta high F1');
save([V1dir '\noTheta_highF1_Ntsr1.mat'],'takeHigh');





a=load([outputDir '\theta_trialAv_temp_noLED.mat']);
[lowF1Trials,highF1Trials]=getHighVsLowF1(a.theta_trialAv_temp,psth,uses,uses_tri,noThetaTrials,0,stimWindow);
firststep=find(noThetaTrials'==0 & ismember(psth.unitLED{1},ledOffVal));

a=load([outputDir '\theta_noLED.mat']);
temp=addAll(a.theta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeLow'); title('Theta low F1');
save([outputDir '\theta_lowF1_dLGN.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials)),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeHigh'); title('Theta high F1');
save([outputDir '\theta_highF1_dLGN.mat'],'takeHigh');

% V1
a=load([V1dir '\theta_noLED.mat']);
temp=addAll(a.theta.allS);

takeLow=reshape(nanmean(temp(:,:,ismember(firststep,lowF1Trials(20:end))),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeLow'); title('V1 units -- Theta low F1');
save([V1dir '\theta_lowF1_Ntsr1.mat'],'takeLow');

takeHigh=reshape(nanmean(temp(:,:,ismember(firststep,highF1Trials(20:end))),3),size(temp,1),size(temp,2));
figure(); imagesc(a.theta.allS.t,a.theta.allS.f,takeHigh'); title('V1 units -- Theta high F1');
save([V1dir '\theta_highF1_Ntsr1.mat'],'takeHigh');


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

function [lowF1Trials,highF1Trials]=getHighVsLowF1(specs_temp,psth,uses,uses_tri,noThetaTrials,isNoTheta,stimWindow)

[lowF1Trials,highF1Trials]=findTrialsWithHighOrLowF1(specs_temp,psth,uses,uses_tri,noThetaTrials,isNoTheta,stimWindow);

end
