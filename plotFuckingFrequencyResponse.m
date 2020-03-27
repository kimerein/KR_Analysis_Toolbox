function [out1,out2]=plotFuckingFrequencyResponse(freqs,noTheta_noLED,noTheta_LED,thresh,allNormVals,allNormVals_LED)

% use thresh 0.1e4 for both no theta and theta

% USE FOR NO THETA
takeTheseUnits=nansum(noTheta_noLED,2)>thresh & nansum(noTheta_LED,2)>thresh;
% takeTheseUnits=nansum(noTheta_noLED,2)>thresh;
normVals_noLED=allNormVals;
% normVals_noLED=allNormVals+noTheta_noLED(:,10)';
% normVals_noLED=allNormVals_LED;
% normVals_noLED=(allNormVals_LED+allNormVals+noTheta_LED(:,10)')/2;
% normVals_LED=(allNormVals_LED+allNormVals+noTheta_LED(:,10)')/2;
normVals_LED=allNormVals_LED;
% normVals_LED=allNormVals_LED+noTheta_LED(:,10)';

% % USE FOR THETA
% % takeTheseUnits=nansum(noTheta_noLED,2)>thresh & nansum(noTheta_LED,2)>thresh;
% takeTheseUnits=nansum(noTheta_noLED,2)>thresh;
% % normVals_noLED=allNormVals;
% % normVals_noLED=allNormVals+noTheta_noLED(:,10)';
% % normVals_noLED=allNormVals_LED;
% normVals_noLED=(allNormVals_LED+allNormVals+noTheta_LED(:,10)')/2;
% normVals_LED=(allNormVals_LED+allNormVals+noTheta_LED(:,10)')/2;
% % normVals_LED=allNormVals_LED;
% % normVals_LED=allNormVals_LED+noTheta_LED(:,10)';

disp('Using this many units');
disp(nansum(takeTheseUnits));

doSmooth=false;
ds=3;
if doSmooth==true
    for i=1:size(noTheta_noLED,1)
        noTheta_noLED(i,:)=smooth(noTheta_noLED(i,:),ds);
        noTheta_LED(i,:)=smooth(noTheta_LED(i,:),ds);
    end
end

figure(); 
out1=noTheta_noLED(takeTheseUnits,:)./repmat(normVals_noLED(takeTheseUnits)',1,15);
plot(freqs,nanmean(noTheta_noLED(takeTheseUnits,:)./repmat(normVals_noLED(takeTheseUnits)',1,15),1),'Color','k'); 
hold on; 
out2=noTheta_LED(takeTheseUnits,:)./repmat(normVals_LED(takeTheseUnits)',1,15);
plot(freqs,nanmean(noTheta_LED(takeTheseUnits,:)./repmat(normVals_LED(takeTheseUnits)',1,15),1),'Color','b');

for i=1:15
    line([freqs(i) freqs(i)],[nanmean(noTheta_noLED(takeTheseUnits,i)./repmat(normVals_noLED(takeTheseUnits)',1,1),1)-nanstd(noTheta_noLED(takeTheseUnits,i)./repmat(normVals_noLED(takeTheseUnits)',1,1),[],1)./sqrt(85) nanmean(noTheta_noLED(takeTheseUnits,i)./repmat(normVals_noLED(takeTheseUnits)',1,1),1)+nanstd(noTheta_noLED(takeTheseUnits,i)./repmat(normVals_noLED(takeTheseUnits)',1,1),[],1)./sqrt(85)],'Color','k');
end


for i=1:15
    line([freqs(i) freqs(i)],[nanmean(noTheta_LED(takeTheseUnits,i)./repmat(normVals_LED(takeTheseUnits)',1,1),1)-nanstd(noTheta_LED(takeTheseUnits,i)./repmat(normVals_LED(takeTheseUnits)',1,1),[],1)./sqrt(85) nanmean(noTheta_LED(takeTheseUnits,i)./repmat(normVals_LED(takeTheseUnits)',1,1),1)+nanstd(noTheta_LED(takeTheseUnits,i)./repmat(normVals_LED(takeTheseUnits)',1,1),[],1)./sqrt(85)],'Color','b');
end

return

mult_noTheta=visResponses_ntsr1_noTheta_noLED'.*test_noTheta_alphaResponse;
guess_noTheta_LED_fund=noTheta_noLED_fund-repmat(nanmean(mult_noTheta,2)',size(noTheta_noLED_fund,1),1)./(1.4924e+06);
figure(); plot(freqs,nanmean(noTheta_noLED_fund,1),'Color','k'); hold on; plot(freqs,nanmean(guess_noTheta_LED_fund,1),'Color','b');

figure(); plot(freqs,nanmean([diff_just_led; diff_just_gtacr2*10],1));



update_following=nanmean(visResponses_ntsr1_theta_noLED,1);
update_following(ismember(freqs,[60]))=update_following(ismember(freqs,[50]));
update_following=update_following+2000;
update_following=smooth(update_following,3);
update_following=update_following';
% set n.s. F1 to 0, n.s. for no-theta are 6,8,10,12,16,18,20,30,40,50,60 Hz
% update_following(ismember(freqs,[6,8,10,12,16,18,20,30,40,50,60]))=0;
% update_following(ismember(freqs,[16,18,20,30,40,50,60]))=0;
figure(); 
plot(freqs,update_following);
update_alpha=nanmean(theta_alphaResponse,2)'-nanmean(theta_alphaResponse(3,:),2)+2000; % range of alpha diff from baseline is 21% of baseline spont alpha, so add appropriate
% update_alpha=nanmean(theta_alphaResponse,2)'; % range of alpha diff from baseline is 21% of baseline spont alpha, so add appropriate
% set n.s. alpha to baseline 17249, n.s. for no-theta are
% 2,4,6,8,10,12,16,18,20,30,40,50,60
% update_alpha(ismember(freqs,[2,4,6,8,10,12,16,18,20,30,40,50,60]))=17249;
% update_alpha(ismember(freqs,[16,18,20,30,40,50,60]))=17249;
% update_alpha(ismember(freqs,[20,30]))=17249;
% update_alpha=smooth(update_alpha,12);
% update_alpha=smooth(update_alpha,1);
% update_alpha(ismember(freqs,[6 8 10 12 14 16 18 20]))=(9500/30).*[6 8 10 12 14 16 18 20];
update_alpha(ismember(freqs,[1 2 4 6 8 10 12 14 16 18 20]))=(9500/30).*[1 2 4 6 8 10 12 14 16 18 20];
update_alpha(ismember(freqs,[30 40 50 60]))=-(9500/30).*[30 40 50 60]+2*(9500/30)+10*(9500/30)+3800+11400;
% update_alpha(ismember(freqs,[60]))=0;
% update_alpha=update_alpha';
figure(); 
plot(freqs,update_alpha);
figure();
plot(freqs,smooth(-update_following.*update_alpha,5),'Color','r');
diff_both=nanmean(all_theta_LED,1)-nanmean(all_theta,1);
diff_both=smooth(diff_both,4);
diff_both=diff_both';
figure(); plot(freqs,diff_both,'Color','k');




update_following=nanmean(visResponses_ntsr1_noTheta_noLED,1);
update_following(ismember(freqs,[14,16,18,20,30,40,50,60]))=interp1([1 2.4],[update_following(ismember(freqs,[14])) update_following(ismember(freqs,[60]))],[1 1.2 1.4 1.6 1.8 2 2.2 2.4]);
update_following=smooth(update_following,3);
update_following=update_following';
% set n.s. F1 to 0, n.s. for no-theta are 6,8,10,12,16,18,20,30,40,50,60 Hz
% update_following(ismember(freqs,[6,8,10,12,16,18,20,30,40,50,60]))=0;
% update_following(ismember(freqs,[16,18,20,30,40,50,60]))=0;
figure(); 
plot(freqs,update_following);
update_alpha=nanmean(noTheta_alphaResponse,2)'+17249; % range of alpha diff from baseline is 21% of baseline spont alpha, so add appropriate
% set n.s. alpha to baseline 17249, n.s. for no-theta are
% 2,4,6,8,10,12,16,18,20,30,40,50,60
% update_alpha(ismember(freqs,[2,4,6,8,10,12,16,18,20,30,40,50,60]))=17249;
% update_alpha(ismember(freqs,[16,18,20,30,40,50,60]))=17249;
% update_alpha(ismember(freqs,[20,30]))=17249;
update_alpha=smooth(update_alpha,4);
update_alpha=update_alpha';
figure(); 
plot(freqs,update_alpha);
figure();
plot(freqs,-update_following.*update_alpha,'Color','r');
% diff_both=nanmean(all_noTheta_LED,1)-nanmean(all_noTheta,1);
diff_just_led=nanmean(noTheta_LED_fund,1)-nanmean(noTheta_noLED_fund,1);
diff_just_led(12:15)=diff_just_led(11);
diff_just_gtacr2=nanmean(noTheta_LED_outForFig,1)-nanmean(noTheta_noLED_outForFig,1);
diff_both=diff_just_gtacr2+0.5*diff_just_led;
diff_both=smooth(diff_both,3);
diff_both=diff_both';
figure(); plot(freqs,diff_both,'Color','k');

% c=-nanmean(visResponses_ntsr1_noTheta_noLED(:,8),1);
c=-0.5*nanmean(visResponses_ntsr1_noTheta_noLED(:,8),1);
% c=-2e4;
% c=-nanmax(nanmean(visResponses_ntsr1_noTheta_noLED,1));
% zeroed_visResponses_ntsr1_noTheta_noLED=nanmean(visResponses_ntsr1_noTheta_noLED,1)+c;
% zeroed_visResponses_ntsr1_noTheta_noLED(9:end)=0;
% c=0;
% update_following=zeroed_visResponses_ntsr1_noTheta_noLED;
update_following=nanmean(visResponses_ntsr1_noTheta_noLED,1)+c;
% update_following(12:end)=0;
figure(); 
plot(freqs,update_following);
% plot(freqs,sign(update_following).*sqrt(abs(update_following)));
figure(); 
update_alpha=nanmean(noTheta_alphaResponse,2)';
% update_alpha=nanmean(test_noTheta_alphaResponse,2)';
update_alpha=update_alpha-nanmin(update_alpha)+0.01e4-2000;
% plot(freqs,sign(update_alpha).*sqrt(abs(update_alpha)));
plot(freqs,update_alpha);
figure();
plot(freqs,-update_following.*update_alpha,'Color','r');
figure(); plot(freqs,diff_both,'Color','k');
figure();
plot(freqs,nanmean(all_noTheta,1));
% plot(freqs,-(sign(update_following).*sqrt(abs(update_following))).*(sign(update_alpha).*sqrt(abs(update_alpha))),'Color','r');
% amp_diff_both=sqrt(nanmean(all_noTheta_LED,1))-sqrt(nanmean(all_noTheta,1));
% figure(); plot(freqs,amp_diff_both,'Color','k');
% figure(); plot(freqs,diff_both,'Color','k');
% feedbackSignal=-(sign(update_following).*sqrt(abs(update_following))).*(sign(update_alpha).*sqrt(abs(update_alpha)));
feedbackSignal=-update_following.*update_alpha+3*10^7;
% ifMultiply=feedbackSignal.*sqrt(nanmean(all_noTheta,1))-sqrt(nanmean(all_noTheta,1));
ifMultiply=(feedbackSignal.*nanmean(all_noTheta,1))./nansum(feedbackSignal.*nanmean(all_noTheta,1))-nanmean(all_noTheta,1)./nansum(nanmean(all_noTheta,1));
figure(); plot(freqs,ifMultiply,'Color','g');
% figure(); plot(freqs,feedbackSignal.*sqrt(nanmean(all_noTheta,1)),'Color','m');
figure(); plot(freqs,feedbackSignal.*nanmean(all_noTheta,1),'Color','m');
% figure(); plot(freqs,feedbackSignal.*nanmean(all_noTheta,1).*nanmean(all_noTheta,1),'Color','m');
% figure(); plot(freqs,feedbackSignal.*nanmean(all_noTheta,1).*nanmean(all_noTheta,1).*nanmean(all_noTheta,1),'Color','m');
