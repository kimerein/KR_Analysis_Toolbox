function [newF1psth_noThetaNoLED,newF1psth_thetaNoLED]=checkForChangeWithDifferentThresh(currdir,thetaThresh,thetaDiff,whichStim,psth,newThetaThresh,F1range)

a=load([currdir '\noTheta_trialAv_temp_noLED.mat']);
noTheta_trialAv_temp=a.noTheta_trialAv_temp;
a=load([currdir '\theta_trialAv_temp_noLED.mat']);
theta_trialAv_temp=a.theta_trialAv_temp;

f=noTheta_trialAv_temp(1).allS.f;
t=noTheta_trialAv_temp(1).allS.t;
l=psth.unitLED{1};
s=psth.unitStimcond{1};

% these were no theta no led
fi=find(nanmean(thetaDiff,2)<thetaThresh & l'==0 & s'==whichStim);
fiNew=find(nanmean(thetaDiff,2)<newThetaThresh & l'==0 & s'==whichStim);
takeThese=ismember(fi,fiNew);
newF1psth_noThetaNoLED=nan(length(noTheta_trialAv_temp(1).allS.S),length(t));
for i=1:length(noTheta_trialAv_temp(1).allS.S)
    spectemp=noTheta_trialAv_temp(1).allS.S{i};
    F1avAcrossNewTrials=nanmean(nanmean(spectemp(:,f>=F1range(1) & f<=F1range(2),takeThese),3),2);
    newF1psth_noThetaNoLED(i,:)=F1avAcrossNewTrials';
end
temp=newF1psth_noThetaNoLED;
evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t<3.75 & t>3),2);
noThetaNoLED_evoked=evoked;

% these were theta no led
fi=find(nanmean(thetaDiff,2)>=thetaThresh & l'==0 & s'==whichStim);
fiNew=find(nanmean(thetaDiff,2)>=newThetaThresh & l'==0 & s'==whichStim);
takeThese=ismember(fi,fiNew);
%takeThese=ones(1,length(fi));
newF1psth_thetaNoLED=nan(length(noTheta_trialAv_temp(1).allS.S),length(t));
for i=1:length(noTheta_trialAv_temp(1).allS.S)
    spectemp=theta_trialAv_temp(1).allS.S{i};
    F1avAcrossNewTrials=nanmean(nanmean(spectemp(:,f>=F1range(1) & f<=F1range(2),takeThese),3),2);
    newF1psth_thetaNoLED(i,:)=F1avAcrossNewTrials';
end
temp=newF1psth_thetaNoLED;
evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t<3.75 & t>3),2);
thetaNoLED_evoked=evoked;

p=signrank(thetaNoLED_evoked-noThetaNoLED_evoked);
disp('no theta no LED vs theta no LED signrank p: ');
disp(p);

figure(); ecdf(thetaNoLED_evoked-noThetaNoLED_evoked,'Bounds','on');

end