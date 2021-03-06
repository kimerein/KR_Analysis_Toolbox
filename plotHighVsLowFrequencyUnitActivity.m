function plotHighVsLowFrequencyUnitActivity(psth,plotThisUnit,noLEDvals,LEDvals,outputDir,noThetaTrials,takeTri)

ds=5;

if ~isempty(takeTri)
    for i=1:length(psth.psths)
        temp=psth.psths{i};
        temp=temp(takeTri,:);
        psth.psths{i}=temp;
        
        temp=psth.unitTrials{i};
        temp=temp(takeTri);
        psth.unitTrials{i}=temp;
        
        temp=psth.unitStimcond{i};
        temp=temp(takeTri);
        psth.unitStimcond{i}=temp;
        
        temp=psth.unitLED{i};
        temp=temp(takeTri);
        psth.unitLED{i}=temp;
    end
    noThetaTrials=noThetaTrials(takeTri);
end

temp=psth.psths{plotThisUnit}; 
figure(); 
plot(downSampAv(psth.t,ds),smooth(downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(noLEDvals)) & noThetaTrials'==1,:),1),ds),3),'Color','k'); 
hold on; 
plot(downSampAv(psth.t,ds),smooth(downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(LEDvals)) & noThetaTrials'==1,:),1),ds),3),'Color','c');
title('No theta');
legend({'led off','led on'});

temp=psth.psths{plotThisUnit}; 
figure(); 
plot(downSampAv(psth.t,ds),smooth(downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(noLEDvals)) & noThetaTrials'==0,:),1),ds),3),'Color','k'); 
hold on; 
plot(downSampAv(psth.t,ds),smooth(downSampAv(nanmean(temp(ismember(single(psth.unitLED{plotThisUnit}),single(LEDvals)) & noThetaTrials'==0,:),1),ds),3),'Color','c');
title('Theta');
legend({'led off','led on'});

return

plotDir=[outputDir '\'];

a=load([outputDir '\noTheta_LED.mat']);
times=a.noTheta.allS.t;

a=load([plotDir 'noTheta_noLED.mat']);
F1_noTheta_noLED=a.noTheta.F1amp;
LFa_noTheta_noLED=a.noTheta.LFa;
HFa_noTheta_noLED=a.noTheta.HFa;
a=load([plotDir 'noTheta_LED.mat']);
F1_noTheta_LED=a.noTheta.F1amp;
a=load([plotDir 'theta_noLED.mat']);
F1_theta_noLED=a.theta.F1amp;
LFa_theta_noLED=a.theta.LFa;
HFa_theta_noLED=a.theta.HFa;
a=load([plotDir 'theta_LED.mat']);
F1_theta_LED=a.theta.F1amp;

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

% figure();
% plot(times,LFa_noTheta_noLED(plotThisUnit,:),'Color','b');
% hold on;
% plot(times,HFa_noTheta_noLED(plotThisUnit,:),'Color','r');
% title('No theta LFa vs HFa');
% legend({'low frequency alpha','high frequency alpha'});
% 
% figure();
% plot(times,HFa_noTheta_noLED(plotThisUnit,:)./LFa_noTheta_noLED(plotThisUnit,:),'Color','g');
% title('No theta ratio HFa to LFa');
% ylim([0 3]);

% figure();
% plot(times,LFa_theta_noLED(plotThisUnit,:),'Color','b');
% hold on;
% plot(times,HFa_theta_noLED(plotThisUnit,:),'Color','r');
% title('Theta LFa vs HFa');
% legend({'low frequency alpha','high frequency alpha'});
% 
% figure();
% plot(times,HFa_theta_noLED(plotThisUnit,:)./LFa_theta_noLED(plotThisUnit,:),'Color','g');
% title('Theta ratio HFa to LFa');
% ylim([0 3]);