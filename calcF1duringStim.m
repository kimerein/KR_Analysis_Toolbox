function summary=calcF1duringStim(output)

base=[0 4];
baseSubtract=1;
stimWindow=[4 6.5];

times=linspace(1,13.5,size(output.noTheta_trialAv_noLED.F1amp,2));

if baseSubtract==0
    noTheta_trialAv_noLED=nanmean(output.noTheta_trialAv_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    theta_trialAv_noLED=nanmean(output.theta_trialAv_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    noTheta_trialAv_LED=nanmean(output.noTheta_trialAv_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    theta_trialAv_LED=nanmean(output.theta_trialAv_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    noTheta_noLED=nanmean(output.noTheta_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    theta_noLED=nanmean(output.theta_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    noTheta_LED=nanmean(output.noTheta_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
    theta_LED=nanmean(output.theta_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2);
else
    noTheta_trialAv_noLED=nanmean(output.noTheta_trialAv_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.noTheta_trialAv_noLED.F1amp(:,times>=base(1) & times<=base(2)),2);
    theta_trialAv_noLED=nanmean(output.theta_trialAv_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.theta_trialAv_noLED.F1amp(:,times>=base(1) & times<=base(2)),2);
    noTheta_trialAv_LED=nanmean(output.noTheta_trialAv_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.noTheta_trialAv_LED.F1amp(:,times>=base(1) & times<=base(2)),2);
    theta_trialAv_LED=nanmean(output.theta_trialAv_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.theta_trialAv_LED.F1amp(:,times>=base(1) & times<=base(2)),2);
    noTheta_noLED=nanmean(output.noTheta_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.noTheta_noLED.F1amp(:,times>=base(1) & times<=base(2)),2);
    theta_noLED=nanmean(output.theta_noLED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.theta_noLED.F1amp(:,times>=base(1) & times<=base(2)),2);
    noTheta_LED=nanmean(output.noTheta_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.noTheta_LED.F1amp(:,times>=base(1) & times<=base(2)),2);
    theta_LED=nanmean(output.theta_LED.F1amp(:,times>=stimWindow(1) & times<=stimWindow(2)),2)-nanmean(output.theta_LED.F1amp(:,times>=base(1) & times<=base(2)),2);
end

figure();
% scatter(noTheta_trialAv_noLED,theta_trialAv_noLED);
% scatter(noTheta_trialAv_noLED,noTheta_trialAv_LED);
scatter(theta_trialAv_noLED,theta_trialAv_LED);
summary.noTheta_trialAv_noLED=noTheta_trialAv_noLED;
summary.theta_trialAv_noLED=theta_trialAv_noLED;
summary.noTheta_trialAv_LED=noTheta_trialAv_LED;
summary.theta_trialAv_LED=theta_trialAv_LED;