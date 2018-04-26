function arrangeF1byFreq(filedir)

a=load([filedir 'noTheta_trialAv_temp_noLED.mat']);
noTheta_noLED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_noLED.mat']);
theta_noLED=a.theta_trialAv_temp;
a=load([filedir 'noTheta_trialAv_temp_LED.mat']);
noTheta_LED=a.noTheta_trialAv_temp;
a=load([filedir 'theta_trialAv_temp_LED.mat']);
theta_LED=a.theta_trialAv_temp;

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
for i=1:length(freqs)
    noTheta_trialAv_noLED=noTheta_noLED(i);
    theta_trialAv_noLED=theta_noLED(i);
    noTheta_trialAv_LED=noTheta_LED(i);
    theta_trialAv_LED=theta_LED(i);
    currdir=[filedir '\Hz' num2str(freqs(i))];
    mkdir(currdir);
    noTheta_trialAv=noTheta_noLED(i);
    theta_trialAv=theta_noLED(i);
    save([currdir '\noTheta_trialAv_noLED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_noLED.mat'],'theta_trialAv');
    noTheta_trialAv=noTheta_LED(i);
    theta_trialAv=theta_LED(i);
    save([currdir '\noTheta_trialAv_LED.mat'],'noTheta_trialAv');
    save([currdir '\theta_trialAv_LED.mat'],'theta_trialAv');
end