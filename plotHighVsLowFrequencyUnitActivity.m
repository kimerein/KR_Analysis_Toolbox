function plotHighVsLowFrequencyUnitActivity(psth,plotThisUnit,noLEDvals,LEDvals,outputDir,noThetaTrials)

temp=psth.psths{plotThisUnit}; 
figure(); 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(noLEDvals)) & noThetaTrials'==1,:),1),10),'Color','k'); 
hold on; 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(LEDvals)) & noThetaTrials'==1,:),1),10),'Color','c');
title('No theta');
legend({'led off','led on'});

temp=psth.psths{plotThisUnit}; 
figure(); 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(noLEDvals)) & noThetaTrials'==0,:),1),10),'Color','k'); 
hold on; 
plot(downSampAv(psth.t,10),downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(LEDvals)) & noThetaTrials'==0,:),1),10),'Color','c');
title('Theta');
legend({'led off','led on'});

plotDir=[outputDir '\'];

a=load([outputDir '\noTheta_LED.mat']);
times=a.noTheta.allS.t;

a=load([plotDir 'noTheta_trialAv_noLED.mat']);
F1_noTheta_noLED=a.noTheta_trialAv.F1amp;
LFa_noTheta_noLED=a.noTheta_trialAv.LFa;
HFa_noTheta_noLED=a.noTheta_trialAv.HFa;
a=load([plotDir 'noTheta_trialAv_LED.mat']);
F1_noTheta_LED=a.noTheta_trialAv.F1amp;
a=load([plotDir 'theta_trialAv_noLED.mat']);
F1_theta_noLED=a.theta_trialAv.F1amp;
LFa_theta_noLED=a.theta_trialAv.LFa;
HFa_theta_noLED=a.theta_trialAv.HFa;
a=load([plotDir 'theta_trialAv_LED.mat']);
F1_theta_LED=a.theta_trialAv.F1amp;

figure();
plot(times,F1_noTheta_noLED(plotThisUnit,:),'Color','k');
hold on;
plot(times,F1_noTheta_LED(plotThisUnit,:),'Color','c');
title('No Theta');

figure();
plot(times,F1_theta_noLED(plotThisUnit,:),'Color','k');
hold on;
plot(times,F1_theta_LED(plotThisUnit,:),'Color','c');
title('Theta');

figure();
plot(times,LFa_noTheta_noLED(plotThisUnit,:),'Color','b');
hold on;
plot(times,HFa_noTheta_noLED(plotThisUnit,:),'Color','r');
title('No theta LFa vs HFa');
legend({'low frequency alpha','high frequency alpha'});

figure();
plot(times,LFa_theta_noLED(plotThisUnit,:),'Color','b');
hold on;
plot(times,HFa_theta_noLED(plotThisUnit,:),'Color','r');
title('Theta LFa vs HFa');
legend({'low frequency alpha','high frequency alpha'});