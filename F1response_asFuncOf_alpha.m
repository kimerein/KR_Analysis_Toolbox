function output=F1response_asFuncOf_alpha(datadir,freqs,normToCon,removeOutliers)

doSmooth=false;

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
    
    if length(a.F1sortedbyalpha.all_noTheta_noLED_F1_sortedByAlpha_low')~=size(noTheta_noLED_F1_low,1)
        disp('skipped something');
        continue
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

% for i=1:15
%     temp=noTheta_noLED_F1_low(:,i);
%     temp(temp<0)=nan;
%     noTheta_noLED_F1_low(:,i)=temp;
%     temp=noTheta_noLED_F1_high(:,i);
%     temp(temp<0)=nan;
%     noTheta_noLED_F1_high(:,i)=temp;
%     
%     temp=noTheta_LED_F1_low(:,i);
%     temp(temp<0)=nan;
%     noTheta_LED_F1_low(:,i)=temp;
%     temp=noTheta_LED_F1_high(:,i);
%     temp(temp<0)=nan;
%     noTheta_LED_F1_high(:,i)=temp;
%     
%     temp=theta_noLED_F1_low(:,i);
%     temp(temp<0)=nan;
%     theta_noLED_F1_low(:,i)=temp;
%     temp=theta_noLED_F1_high(:,i);
%     temp(temp<0)=nan;
%     theta_noLED_F1_high(:,i)=temp;
%     
%     temp=theta_LED_F1_low(:,i);
%     temp(temp<0)=nan;
%     theta_LED_F1_low(:,i)=temp;
%     temp=theta_LED_F1_high(:,i);
%     temp(temp<0)=nan;
%     theta_LED_F1_high(:,i)=temp;
% end

if removeOutliers==true
    saway=2.2;
    temp=noTheta_noLED_F1_low(:,5);
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_noLED_F1_low(:,5)=temp;
    
    temp=noTheta_noLED_F1_high(:,5);
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_noLED_F1_high(:,5)=temp;
    
%     saway=2.2;
%     temp=theta_noLED_F1_low(:,1);
%     mes=nanmean(temp,1);
%     ss=nanstd(temp,[],1);
%     maxInRange=mes+saway*ss;
%     minInRange=mes-saway*ss;
%     temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
%     temp(temp<repmat(minInRange,size(temp,1),1))=nan;
%     theta_noLED_F1_low(:,1)=temp;
%     
%     temp=theta_noLED_F1_high(:,1);
%     mes=nanmean(temp,1);
%     ss=nanstd(temp,[],1);
%     maxInRange=mes+saway*ss;
%     minInRange=mes-saway*ss;
%     temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
%     temp(temp<repmat(minInRange,size(temp,1),1))=nan;
%     theta_noLED_F1_high(:,1)=temp;
    
    saway=2.2;
    temp=theta_noLED_F1_low(:,10);
%     mes=nanmean(temp,1);
%     ss=nanstd(temp,[],1);
%     maxInRange=mes+saway*ss;
%     minInRange=mes-saway*ss;
%     temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
%     temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    temp(temp<12)=nan;
    theta_noLED_F1_low(:,10)=temp;
    
    temp=theta_noLED_F1_high(:,10);
%     mes=nanmean(temp,1);
%     ss=nanstd(temp,[],1);
%     maxInRange=mes+saway*ss;
%     minInRange=mes-saway*ss;
%     temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
%     temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    temp(temp<12)=nan;
    theta_noLED_F1_high(:,10)=temp;
end
    
if removeOutliers==true
%     saway=3.4;
%     saway=3.2;
    saway=5.6;
%     saway=2.6;
    
    temp=noTheta_noLED_F1_low;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_noLED_F1_low=temp;
    
    temp=noTheta_noLED_F1_high;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_noLED_F1_high=temp;
    
    temp=theta_noLED_F1_low;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    theta_noLED_F1_low=temp;
    
    temp=theta_noLED_F1_high;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    theta_noLED_F1_high=temp;
    
    
    temp=noTheta_LED_F1_low;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_LED_F1_low=temp;
    
    temp=noTheta_LED_F1_high;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    noTheta_LED_F1_high=temp;
    
    temp=theta_LED_F1_low;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    theta_LED_F1_low=temp;
    
    temp=theta_LED_F1_high;
    mes=nanmean(temp,1);
    ss=nanstd(temp,[],1);
    maxInRange=mes+saway*ss;
    minInRange=mes-saway*ss;
    temp(temp>repmat(maxInRange,size(temp,1),1))=nan;
    temp(temp<repmat(minInRange,size(temp,1),1))=nan;
    theta_LED_F1_high=temp;
end

if doSmooth==true
    for i=1:size(noTheta_noLED_F1_low,1)
        noTheta_noLED_F1_low(i,5:end)=smooth(noTheta_noLED_F1_low(i,5:end),3);
        noTheta_noLED_F1_high(i,5:end)=smooth(noTheta_noLED_F1_high(i,5:end),3);

        noTheta_LED_F1_low(i,5:end)=smooth(noTheta_LED_F1_low(i,5:end),3);
        noTheta_LED_F1_high(i,5:end)=smooth(noTheta_LED_F1_high(i,5:end),3);
        
        theta_noLED_F1_low(i,5:end)=smooth(theta_noLED_F1_low(i,5:end),3);
        theta_noLED_F1_high(i,5:end)=smooth(theta_noLED_F1_high(i,5:end),3);
        
        theta_LED_F1_low(i,5:end)=smooth(theta_LED_F1_low(i,5:end),3);
        theta_LED_F1_high(i,5:end)=smooth(theta_LED_F1_high(i,5:end),3);
    end
end


if normToCon==true
    normFac=nanmean(nanmean(noTheta_noLED_F1_high));
else
    normFac=1;
end
figure();
for i=1:15
    temp=noTheta_noLED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp)./normFac,[],'k');
end
hold on;
for i=1:15
    temp=noTheta_noLED_F1_high(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','b');
    scatter(freqs(i),nanmean(temp)./normFac,[],'r');
end
set(gca,'Xscale','log');
title('No theta LOW ALPHA is Black');

if normToCon==true
    normFac=nanmean(nanmean(theta_noLED_F1_high));
    output.normFac=normFac;
else
    normFac=1;
end
figure();
for i=1:15
    temp=theta_noLED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp)./normFac,[],'k');
end
hold on;
for i=1:15
    temp=theta_noLED_F1_high(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','b');
    scatter(freqs(i),nanmean(temp)./normFac,[],'r');
end
set(gca,'Xscale','log');
title('Theta LOW ALPHA is Black');

if normToCon==true
    normFac=nanmean(nanmean(noTheta_noLED_F1_high));
else
    normFac=1;
end
figure();
for i=1:15
    temp=noTheta_noLED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp)./normFac,[],'k');
end
hold on;
for i=1:15
%     temp=(noTheta_LED_F1_low(:,i)+noTheta_LED_F1_high(:,i))/2;
    temp=noTheta_LED_F1_low(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','b');
    scatter(freqs(i),nanmean(temp)./normFac,[],'r');
end
set(gca,'Xscale','log');
title('No theta LED is Red');

if normToCon==true
    normFac=nanmean(nanmean(theta_noLED_F1_high));
else
    normFac=1;
end
figure();
for i=1:15
    temp=(theta_noLED_F1_low(:,i)+theta_noLED_F1_high(:,i))/2;
%     temp=theta_noLED_F1_high(:,i);
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','k');
    hold on;
    scatter(freqs(i),nanmean(temp)./normFac,[],'k');
end
hold on;
for i=1:15
    temp=(theta_LED_F1_low(:,i)+theta_LED_F1_high(:,i))/2;
    line([freqs(i) freqs(i)],[nanmean(temp)-nanstd(temp,[],1)./sqrt(length(temp)) nanmean(temp)+nanstd(temp,[],1)./sqrt(length(temp))]./normFac,'Color','b');
    scatter(freqs(i),nanmean(temp)./normFac,[],'r');
end
set(gca,'Xscale','log');
title('Theta LED is Red');

output.set1=theta_noLED_F1_low;
output.set2=theta_noLED_F1_high;