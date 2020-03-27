function [noTheta_noLED_out,noTheta_LED_out,theta_noLED_out,theta_LED_out]=getF1distributions(dd,takeTrials_noTheta,takeTrials_noTheta_LED,takeTrials_theta,takeTrials_theta_LED,F1range,alphaRange,timeWindow,baseWindow,suppressOutput,baseTrials1,baseTrials2,baseTrials3,baseTrials4,otherWindow)

baseIsSameTimeWindow=true;
getAlphaFromOther=true;

% no-theta

a=load([dd '\noTheta_noLED.mat']);
noTheta=a.noTheta;

[~,F1s]=plotF1distributionsZscored(noTheta,takeTrials_noTheta,F1range,timeWindow,true);
if getAlphaFromOther==false
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_noTheta,alphaRange,timeWindow,true);
else
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_noTheta,alphaRange,otherWindow,true);
end
if baseIsSameTimeWindow==false
    [~,F1s_base]=plotF1distributionsZscored(noTheta,takeTrials_noTheta,F1range,baseWindow,true);
else
    if isempty(baseTrials1)
        temp=1:size(noTheta.allS.S{1},3);
        [~,F1s_base]=plotF1distributionsZscored(noTheta,temp(~ismember(temp,takeTrials_noTheta)),F1range,timeWindow,true);
    else
        [~,F1s_base]=plotF1distributionsZscored(noTheta,baseTrials1,F1range,timeWindow,true);
    end
end
[~,alphas_base]=plotF1distributionsZscored(noTheta,takeTrials_noTheta,alphaRange,baseWindow,true);
if baseIsSameTimeWindow==false
    F1s=F1s-F1s_base;
else
    F1s=F1s-repmat(nanmean(F1s_base,2),1,size(F1s,2));
end
% alphas=alphas-alphas_base;

isResponsive=nanmean(F1s,2)>nanmedian(nanmean(F1s,2),1);
temp=F1s;
temp_alpha=alphas;
% temp=temp(isResponsive,:);
% temp_alpha=temp_alpha(isResponsive,:);
temp=temp(1:end); 
temp_alpha=temp_alpha(1:end);

if suppressOutput==false
    [n,x]=hist(temp,50);
    figure();
    plot(x,n,'Color','k');
    hold on;
    [n,x]=hist(temp(temp_alpha<1));
    plot(x,n,'Color','g');
    [n,x]=hist(temp(temp_alpha>1));
    plot(x,n,'Color','r');
    title('No theta No LED');
end

noTheta_noLED_out.F1s=F1s;
noTheta_noLED_out.alphas=alphas;

% no-theta LED
a=load([dd '\noTheta_LED.mat']);
noTheta=a.noTheta;

[~,F1s]=plotF1distributionsZscored(noTheta,takeTrials_noTheta_LED,F1range,timeWindow,true);
if getAlphaFromOther==false    
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_noTheta_LED,alphaRange,timeWindow,true);
else
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_noTheta_LED,alphaRange,otherWindow,true);
end
if baseIsSameTimeWindow==false
    [~,F1s_base]=plotF1distributionsZscored(noTheta,takeTrials_noTheta_LED,F1range,baseWindow,true);
else
    if isempty(baseTrials2)
        temp=1:size(noTheta.allS.S{1},3);
        [~,F1s_base]=plotF1distributionsZscored(noTheta,temp(~ismember(temp,takeTrials_noTheta_LED)),F1range,timeWindow,true);
    else
        [~,F1s_base]=plotF1distributionsZscored(noTheta,baseTrials2,F1range,timeWindow,true);
    end
end
[~,alphas_base]=plotF1distributionsZscored(noTheta,takeTrials_noTheta_LED,alphaRange,baseWindow,true);
if baseIsSameTimeWindow==false
    F1s=F1s-F1s_base;
else
    F1s=F1s-repmat(nanmean(F1s_base,2),1,size(F1s,2));
end
% alphas=alphas-alphas_base;

isResponsive=nanmean(F1s,2)>nanmedian(nanmean(F1s,2),1);
temp=F1s;
temp_alpha=alphas;
% temp=temp(isResponsive,:);
% temp_alpha=temp_alpha(isResponsive,:);
temp=temp(1:end); 
temp_alpha=temp_alpha(1:end);

if suppressOutput==false
    [n,x]=hist(temp,50);
    figure();
    plot(x,n,'Color','k');
    hold on;
    [n,x]=hist(temp(temp_alpha<1));
    plot(x,n,'Color','g');
    [n,x]=hist(temp(temp_alpha>1));
    plot(x,n,'Color','r');
    title('No theta LED');
end

noTheta_LED_out.F1s=F1s;
noTheta_LED_out.alphas=alphas;


% theta

a=load([dd '\theta_noLED.mat']);
theta=a.theta;

[~,F1s]=plotF1distributionsZscored(theta,takeTrials_theta,F1range,timeWindow,true);
if getAlphaFromOther==false    
    [~,alphas]=plotF1distributionsZscored(theta,takeTrials_theta,alphaRange,timeWindow,true);
else
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_theta,alphaRange,otherWindow,true);
end
if baseIsSameTimeWindow==false
    [~,F1s_base]=plotF1distributionsZscored(theta,takeTrials_theta,F1range,baseWindow,true);
else
    if isempty(baseTrials3)
        temp=1:size(theta.allS.S{1},3);
        [~,F1s_base]=plotF1distributionsZscored(theta,temp(~ismember(temp,takeTrials_theta)),F1range,timeWindow,true);
    else
        [~,F1s_base]=plotF1distributionsZscored(theta,baseTrials3,F1range,timeWindow,true);
    end
end
[~,alphas_base]=plotF1distributionsZscored(theta,takeTrials_theta,alphaRange,baseWindow,true);
if baseIsSameTimeWindow==false
    F1s=F1s-F1s_base;
else
    F1s=F1s-repmat(nanmean(F1s_base,2),1,size(F1s,2));
end
% alphas=alphas-alphas_base;

isResponsive=nanmean(F1s,2)>nanmedian(nanmean(F1s,2),1);
temp=F1s;
temp_alpha=alphas;
% temp=temp(isResponsive,:);
% temp_alpha=temp_alpha(isResponsive,:);
temp=temp(1:end); 
temp_alpha=temp_alpha(1:end);

if suppressOutput==false
    [n,x]=hist(temp,50);
    figure();
    plot(x,n,'Color','k');
    hold on;
    [n,x]=hist(temp(temp_alpha<1));
    plot(x,n,'Color','g');
    [n,x]=hist(temp(temp_alpha>1));
    plot(x,n,'Color','r');
    title('Theta No LED');
end

theta_noLED_out.F1s=F1s;
theta_noLED_out.alphas=alphas;

% theta LED
a=load([dd '\theta_LED.mat']);
theta=a.theta;

[~,F1s]=plotF1distributionsZscored(theta,takeTrials_theta_LED,F1range,timeWindow,true);
if getAlphaFromOther==false    
    [~,alphas]=plotF1distributionsZscored(theta,takeTrials_theta_LED,alphaRange,timeWindow,true);
else
    [~,alphas]=plotF1distributionsZscored(noTheta,takeTrials_theta_LED,alphaRange,otherWindow,true);
end
if baseIsSameTimeWindow==false
    [~,F1s_base]=plotF1distributionsZscored(theta,takeTrials_theta_LED,F1range,baseWindow,true);
else
    if isempty(baseTrials4)
        temp=1:size(theta.allS.S{1},3);
        [~,F1s_base]=plotF1distributionsZscored(theta,temp(~ismember(temp,takeTrials_theta_LED)),F1range,timeWindow,true);
    else
        [~,F1s_base]=plotF1distributionsZscored(theta,baseTrials4,F1range,timeWindow,true);
    end
end
[~,alphas_base]=plotF1distributionsZscored(theta,takeTrials_theta_LED,alphaRange,baseWindow,true);
if baseIsSameTimeWindow==false
    F1s=F1s-F1s_base;
else
    F1s=F1s-repmat(nanmean(F1s_base,2),1,size(F1s,2));
end
% alphas=alphas-alphas_base;

isResponsive=nanmean(F1s,2)>nanmedian(nanmean(F1s,2),1);
temp=F1s;
temp_alpha=alphas;
% temp=temp(isResponsive,:);
% temp_alpha=temp_alpha(isResponsive,:);
temp=temp(1:end); 
temp_alpha=temp_alpha(1:end);

if suppressOutput==false
    [n,x]=hist(temp,50);
    figure();
    plot(x,n,'Color','k');
    hold on;
    [n,x]=hist(temp(temp_alpha<1));
    plot(x,n,'Color','g');
    [n,x]=hist(temp(temp_alpha>1));
    plot(x,n,'Color','r');
    title('Theta LED');
end

theta_LED_out.F1s=F1s;
theta_LED_out.alphas=alphas;