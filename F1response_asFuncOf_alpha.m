function F1response_asFuncOf_alpha(datadir,freqs,norm_noThetaNoLED,norm_noThetaLED,norm_thetaNoLED,norm_thetaLED)

for i=1:length(freqs)
    a=load([datadir '\F1sorted_' num2str(freqs(i)) 'Hz.mat']);
    
    if i==1
        noTheta_noLED_F1_low=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        noTheta_noLED_F1_high=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        noTheta_LED_F1_low=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        noTheta_LED_F1_high=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        theta_noLED_F1_low=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        theta_noLED_F1_high=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        theta_LED_F1_low=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
        theta_LED_F1_high=nan(length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low),length(freqs));
    end
    
    noTheta_noLED_F1_low(:,i)=a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low';
    noTheta_noLED_F1_high(:,i)=a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_high';
    
    noTheta_LED_F1_low(:,i)=a.F1sortedbyalpha.all_noTheta_LED_F1_sortedByAlpha_low';
    noTheta_LED_F1_high(:,i)=a.F1sortedbyalpha.all_noTheta_LED_F1_sortedByAlpha_high';
    
    theta_noLED_F1_low(:,i)=a.F1sortedbyalpha.all_theta_noLED_F1_sortedByAlpha_low';
    theta_noLED_F1_high(:,i)=a.F1sortedbyalpha.all_theta_noLED_F1_sortedByAlpha_high';
    
    theta_LED_F1_low(:,i)=a.F1sortedbyalpha.all_theta_LED_F1_sortedByAlpha_low';
    theta_LED_F1_high(:,i)=a.F1sortedbyalpha.all_theta_LED_F1_sortedByAlpha_high';
end

figure();
for i=1:15
    temp=noTheta_noLED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))],'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp),[],'k');
end
hold on;
for i=1:15
    temp=noTheta_noLED_F1_high(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))],'Color','b');
    scatter(freqs(i),nanmean(temp),[],'r');
end
set(gca,'Xscale','log');
title('No theta LOW ALPHA is Black');

figure();
for i=1:15
    temp=theta_noLED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))],'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp),[],'k');
end
hold on;
for i=1:15
    temp=theta_noLED_F1_high(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))],'Color','b');
    scatter(freqs(i),nanmean(temp),[],'r');
end
set(gca,'Xscale','log');
title('Theta LOW ALPHA is Black');