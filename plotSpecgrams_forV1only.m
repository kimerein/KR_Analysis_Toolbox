function plotSpecgrams_forV1only(datadir,outputDir)

if iscell(datadir)
%     all_noTheta_trialAv_noLED.t=[];
%     all_noTheta_trialAv_noLED.f=[];
%     all_noTheta_trialAv_noLED.low.S=[];
%     all_noTheta_trialAv_noLED.high.S=[];
%     all_noTheta_trialAv_LED.t=[];
%     all_noTheta_trialAv_LED.f=[];
%     all_noTheta_trialAv_LED.low.S=[];
%     all_noTheta_trialAv_LED.high.S=[];
%     all_theta_trialAv_noLED.t=[];
%     all_theta_trialAv_noLED.f=[];
%     all_theta_trialAv_noLED.low.S=[];
%     all_theta_trialAv_noLED.high.S=[];
%     all_theta_trialAv_LED.t=[];
%     all_theta_trialAv_LED.f=[];
%     all_theta_trialAv_LED.low.S=[];
%     all_theta_trialAv_LED.high.S=[];
    all_allTheta_trialAv_noLED.t=[];
    all_allTheta_trialAv_noLED.f=[];
    all_allTheta_trialAv_noLED.low.S=[];
    all_allTheta_trialAv_noLED.high.S=[];
    all_allTheta_trialAv_LED.t=[];
    all_allTheta_trialAv_LED.f=[];
    all_allTheta_trialAv_LED.low.S=[];
    all_allTheta_trialAv_LED.high.S=[];
    for i=1:length(datadir)
        d=datadir{i};
%         a=load([d '\' 'noTheta_trialAv_temp_noLED']);
%         noTheta_trialAv_noLED=a.noTheta_trialAv_temp;
%         [lowF1_noTheta_trialAv_noLED,highF1_noTheta_trialAv_noLED]=sortLowAndHighF1(noTheta_trialAv_noLED,outputDir,'noTheta_trialAv_noLED');
%         save([outputDir '\' 'lowF1_noTheta_trialAv_noLED.mat'],'lowF1_noTheta_trialAv_noLED');
%         save([outputDir '\' 'highF1_noTheta_trialAv_noLED.mat'],'highF1_noTheta_trialAv_noLED');
%         avSpecgram_low=avS(lowF1_noTheta_trialAv_noLED);
%         avSpecgram_high=avS(highF1_noTheta_trialAv_noLED);
%         if i==1
%             all_noTheta_trialAv_noLED.t=avSpecgram_low.t;
%             all_noTheta_trialAv_noLED.f=avSpecgram_low.f;
%         end
%         all_noTheta_trialAv_noLED.low.S=[all_noTheta_trialAv_noLED.low.S avSpecgram_low.S];
%         all_noTheta_trialAv_noLED.high.S=[all_noTheta_trialAv_noLED.high.S avSpecgram_high.S];
        
%         a=load([d '\' 'noTheta_trialAv_temp_LED']);
%         noTheta_trialAv_LED=a.noTheta_trialAv_temp;
%         [lowF1_noTheta_trialAv_LED,highF1_noTheta_trialAv_LED]=sortLowAndHighF1(noTheta_trialAv_LED,outputDir,'noTheta_trialAv_LED');
%         save([outputDir '\' 'lowF1_noTheta_trialAv_LED.mat'],'lowF1_noTheta_trialAv_LED');
%         save([outputDir '\' 'highF1_noTheta_trialAv_LED.mat'],'highF1_noTheta_trialAv_LED');
%         avSpecgram_low=avS(lowF1_noTheta_trialAv_LED);
%         avSpecgram_high=avS(highF1_noTheta_trialAv_LED);
%         if i==1
%             all_noTheta_trialAv_LED.t=avSpecgram_low.t;
%             all_noTheta_trialAv_LED.f=avSpecgram_low.f;
%         end
%         all_noTheta_trialAv_LED.low.S=[all_noTheta_trialAv_LED.low.S avSpecgram_low.S];
%         all_noTheta_trialAv_LED.high.S=[all_noTheta_trialAv_LED.high.S avSpecgram_high.S];
        
%         a=load([d '\' 'theta_trialAv_temp_noLED']);
%         theta_trialAv_noLED=a.theta_trialAv_temp;
%         [lowF1_theta_trialAv_noLED,highF1_theta_trialAv_noLED]=sortLowAndHighF1(theta_trialAv_noLED,outputDir,'theta_trialAv_noLED');
%         save([outputDir '\' 'lowF1_theta_trialAv_noLED.mat'],'lowF1_theta_trialAv_noLED');
%         save([outputDir '\' 'highF1_theta_trialAv_noLED.mat'],'highF1_theta_trialAv_noLED');
%         avSpecgram_low=avS(lowF1_theta_trialAv_noLED);
%         avSpecgram_high=avS(highF1_theta_trialAv_noLED);
%         if i==1
%             all_theta_trialAv_noLED.t=avSpecgram_low.t;
%             all_theta_trialAv_noLED.f=avSpecgram_low.f;
%         end
%         all_theta_trialAv_noLED.low.S=[all_theta_trialAv_noLED.low.S avSpecgram_low.S];
%         all_theta_trialAv_noLED.high.S=[all_theta_trialAv_noLED.high.S avSpecgram_high.S];
        
%         a=load([d '\' 'theta_trialAv_temp_LED']);
%         theta_trialAv_LED=a.theta_trialAv_temp;
%         [lowF1_theta_trialAv_LED,highF1_theta_trialAv_LED]=sortLowAndHighF1(theta_trialAv_LED,outputDir,'theta_trialAv_LED');
%         save([outputDir '\' 'lowF1_theta_trialAv_LED.mat'],'lowF1_theta_trialAv_LED');
%         save([outputDir '\' 'highF1_theta_trialAv_LED.mat'],'highF1_theta_trialAv_LED');
%         avSpecgram_low=avS(lowF1_theta_trialAv_LED);
%         avSpecgram_high=avS(highF1_theta_trialAv_LED);
%         if i==1
%             all_theta_trialAv_LED.t=avSpecgram_low.t;
%             all_theta_trialAv_LED.f=avSpecgram_low.f;
%         end
%         all_theta_trialAv_LED.low.S=[all_theta_trialAv_LED.low.S avSpecgram_low.S];
%         all_theta_trialAv_LED.high.S=[all_theta_trialAv_LED.high.S avSpecgram_high.S];
        
        a=load([d '\' 'alltheta_trialAv_temp_noLED']);
        alltheta_trialAv_noLED=a.alltheta_trialAv_temp;
%         [lowF1_alltheta_trialAv_noLED,highF1_alltheta_trialAv_noLED]=sortLowAndHighF1(alltheta_trialAv_noLED,outputDir,'alltheta_trialAv_noLED');
%         save([outputDir '\' 'lowF1_alltheta_trialAv_noLED.mat'],'lowF1_alltheta_trialAv_noLED');
%         save([outputDir '\' 'highF1_alltheta_trialAv_noLED.mat'],'highF1_alltheta_trialAv_noLED');
        avSpecgram_low=avS(alltheta_trialAv_noLED);
%         avSpecgram_high=avS(highF1_alltheta_trialAv_noLED);
        if i==1
            all_allTheta_trialAv_noLED.t=avSpecgram_low.t;
            all_allTheta_trialAv_noLED.f=avSpecgram_low.f;
        end
        all_allTheta_trialAv_noLED.low.S=[all_allTheta_trialAv_noLED.low.S avSpecgram_low.S];
%         all_allTheta_trialAv_noLED.high.S=[all_allTheta_trialAv_noLED.high.S avSpecgram_high.S];
        
        a=load([d '\' 'alltheta_trialAv_temp_LED']);
        alltheta_trialAv_LED=a.alltheta_trialAv_temp;
%         [lowF1_alltheta_trialAv_LED,highF1_alltheta_trialAv_LED]=sortLowAndHighF1(alltheta_trialAv_LED,outputDir,'alltheta_trialAv_LED');
%         save([outputDir '\' 'lowF1_alltheta_trialAv_LED.mat'],'lowF1_alltheta_trialAv_LED');
%         save([outputDir '\' 'highF1_alltheta_trialAv_LED.mat'],'highF1_alltheta_trialAv_LED');
        avSpecgram_low=avS(alltheta_trialAv_LED);
%         avSpecgram_high=avS(highF1_alltheta_trialAv_LED);
        if i==1
            all_allTheta_trialAv_LED.t=avSpecgram_low.t;
            all_allTheta_trialAv_LED.f=avSpecgram_low.f;
        end
        all_allTheta_trialAv_LED.low.S=[all_allTheta_trialAv_LED.low.S avSpecgram_low.S];
%         all_allTheta_trialAv_LED.high.S=[all_allTheta_trialAv_LED.high.S avSpecgram_high.S];

    end
%     noTheta_trialAv_noLED=all_noTheta_trialAv_noLED;
%     noTheta_trialAv_LED=all_noTheta_trialAv_LED;
%     theta_trialAv_noLED=all_theta_trialAv_noLED;
%     theta_trialAv_LED=all_theta_trialAv_LED;
    alltheta_trialAv_noLED=all_allTheta_trialAv_noLED;
    alltheta_trialAv_LED=all_allTheta_trialAv_LED;
end












% sumS=zeros(size(noTheta_trialAv_noLED.low.S{1}));
% tally=0;
% for i=1:length(noTheta_trialAv_noLED.low.S)
%     if ~isnan(noTheta_trialAv_noLED.low.S{i})
%         sumS=sumS+noTheta_trialAv_noLED.low.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(noTheta_trialAv_noLED.t,noTheta_trialAv_noLED.f(noTheta_trialAv_noLED.f<=50),sumS(:,noTheta_trialAv_noLED.f<=30)');
% title('noTheta_trialAv_noLED LOW F1');
% save([outputDir '\' 'noTheta_trialAv_noLED LOW F1.mat'],'sumS');
% HFa_lowF1=nanmean(sumS(:,noTheta_trialAv_noLED.f>=11 & noTheta_trialAv_noLED.f<=20),2);
% sumS=zeros(size(noTheta_trialAv_noLED.high.S{1}));
% tally=0;
% for i=1:length(noTheta_trialAv_noLED.high.S)
%     if ~isnan(noTheta_trialAv_noLED.high.S{i})
%         sumS=sumS+noTheta_trialAv_noLED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(noTheta_trialAv_noLED.t,noTheta_trialAv_noLED.f(noTheta_trialAv_noLED.f<=50),sumS(:,noTheta_trialAv_noLED.f<=30)');
% title('noTheta_trialAv_noLED HIGH F1');
% save([outputDir '\' 'noTheta_trialAv_noLED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,noTheta_trialAv_noLED.f>=11 & noTheta_trialAv_noLED.f<=20),2);
% figure(); plot(noTheta_trialAv_noLED.t,HFa_lowF1,'Color','k'); hold on; plot(noTheta_trialAv_noLED.t,HFa_highF1,'Color','r');


% sumS=zeros(size(noTheta_trialAv_LED.low.S{1}));
% tally=0;
% for i=1:length(noTheta_trialAv_LED.low.S)
%     if ~isnan(noTheta_trialAv_LED.low.S{i})
%         sumS=sumS+noTheta_trialAv_LED.low.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(noTheta_trialAv_LED.t,noTheta_trialAv_LED.f(noTheta_trialAv_LED.f<=50),sumS(:,noTheta_trialAv_LED.f<=30)');
% title('noTheta_trialAv_LED LOW F1');
% save([outputDir '\' 'noTheta_trialAv_LED LOW F1.mat'],'sumS');
% HFa_lowF1=nanmean(sumS(:,noTheta_trialAv_LED.f>=11 & noTheta_trialAv_LED.f<=20),2);
% sumS=zeros(size(noTheta_trialAv_LED.high.S{1}));
% tally=0;
% for i=1:length(noTheta_trialAv_LED.high.S)
%     if ~isnan(noTheta_trialAv_LED.high.S{i})
%         sumS=sumS+noTheta_trialAv_LED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(noTheta_trialAv_LED.t,noTheta_trialAv_LED.f(noTheta_trialAv_LED.f<=50),sumS(:,noTheta_trialAv_LED.f<=30)');
% title('noTheta_trialAv_LED HIGH F1');
% save([outputDir '\' 'noTheta_trialAv_LED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,noTheta_trialAv_LED.f>=11 & noTheta_trialAv_LED.f<=20),2);
% figure(); plot(noTheta_trialAv_LED.t,HFa_lowF1,'Color','k'); hold on; plot(noTheta_trialAv_LED.t,HFa_highF1,'Color','r');
% 
% sumS=zeros(size(theta_trialAv_noLED.low.S{1}));
% tally=0;
% for i=1:length(theta_trialAv_noLED.low.S)
%     if ~isnan(theta_trialAv_noLED.low.S{i})
%         sumS=sumS+theta_trialAv_noLED.low.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(theta_trialAv_noLED.t,theta_trialAv_noLED.f(theta_trialAv_noLED.f<=50),sumS(:,theta_trialAv_noLED.f<=30)');
% title('theta_trialAv_noLED LOW F1');
% save([outputDir '\' 'theta_trialAv_noLED LOW F1.mat'],'sumS');
% HFa_lowF1=nanmean(sumS(:,theta_trialAv_noLED.f>=11 & theta_trialAv_noLED.f<=20),2);
% sumS=zeros(size(theta_trialAv_noLED.high.S{1}));
% tally=0;
% for i=1:length(theta_trialAv_noLED.high.S)
%     if ~isnan(theta_trialAv_noLED.high.S{i})
%         sumS=sumS+theta_trialAv_noLED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(theta_trialAv_noLED.t,theta_trialAv_noLED.f(theta_trialAv_noLED.f<=50),sumS(:,theta_trialAv_noLED.f<=30)');
% title('theta_trialAv_noLED HIGH F1');
% save([outputDir '\' 'theta_trialAv_noLED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,theta_trialAv_noLED.f>=11 & theta_trialAv_noLED.f<=20),2);
% figure(); plot(theta_trialAv_noLED.t,HFa_lowF1,'Color','k'); hold on; plot(theta_trialAv_noLED.t,HFa_highF1,'Color','r');
% 
% 
% sumS=zeros(size(theta_trialAv_LED.low.S{1}));
% tally=0;
% for i=1:length(theta_trialAv_LED.low.S)
%     if ~isnan(theta_trialAv_LED.low.S{i})
%         sumS=sumS+theta_trialAv_LED.low.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(theta_trialAv_LED.t,theta_trialAv_LED.f(theta_trialAv_LED.f<=50),sumS(:,theta_trialAv_LED.f<=30)');
% title('theta_trialAv_LED LOW F1');
% save([outputDir '\' 'theta_trialAv_LED LOW F1.mat'],'sumS');
% HFa_lowF1=nanmean(sumS(:,theta_trialAv_LED.f>=11 & theta_trialAv_LED.f<=20),2);
% sumS=zeros(size(theta_trialAv_LED.high.S{1}));
% tally=0;
% for i=1:length(theta_trialAv_LED.high.S)
%     if ~isnan(theta_trialAv_LED.high.S{i})
%         sumS=sumS+theta_trialAv_LED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(theta_trialAv_LED.t,theta_trialAv_LED.f(theta_trialAv_LED.f<=50),sumS(:,theta_trialAv_LED.f<=30)');
% title('theta_trialAv_LED HIGH F1');
% save([outputDir '\' 'theta_trialAv_LED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,theta_trialAv_LED.f>=11 & theta_trialAv_LED.f<=20),2);
% figure(); plot(theta_trialAv_LED.t,HFa_lowF1,'Color','k'); hold on; plot(theta_trialAv_LED.t,HFa_highF1,'Color','r');









sumS=zeros(size(alltheta_trialAv_noLED.low.S{1}));
tally=0;
for i=1:length(alltheta_trialAv_noLED.low.S)
    if ~isnan(alltheta_trialAv_noLED.low.S{i})
        sumS=sumS+alltheta_trialAv_noLED.low.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(alltheta_trialAv_noLED.t,alltheta_trialAv_noLED.f(alltheta_trialAv_noLED.f<=50),sumS(:,alltheta_trialAv_noLED.f<=50)');
title('alltheta_trialAv_noLED ALL F1');
save([outputDir '\' 'alltheta_trialAv_noLED ALL F1.mat'],'sumS');
save([outputDir '\' 'UbyUalltheta_trialAv_noLED.mat'],'alltheta_trialAv_noLED');
% HFa_lowF1=nanmean(sumS(:,alltheta_trialAv_noLED.f>=11 & alltheta_trialAv_noLED.f<=20),2);
% sumS=zeros(size(alltheta_trialAv_noLED.high.S{1}));
% tally=0;
% for i=1:length(alltheta_trialAv_noLED.high.S)
%     if ~isnan(alltheta_trialAv_noLED.high.S{i})
%         sumS=sumS+alltheta_trialAv_noLED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(alltheta_trialAv_noLED.t,alltheta_trialAv_noLED.f(alltheta_trialAv_noLED.f<=50),sumS(:,alltheta_trialAv_noLED.f<=30)');
% title('alltheta_trialAv_noLED HIGH F1');
% save([outputDir '\' 'alltheta_trialAv_noLED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,alltheta_trialAv_noLED.f>=11 & alltheta_trialAv_noLED.f<=20),2);
% figure(); plot(alltheta_trialAv_noLED.t,HFa_lowF1,'Color','k'); hold on; plot(alltheta_trialAv_noLED.t,HFa_highF1,'Color','r');


sumS=zeros(size(alltheta_trialAv_LED.low.S{1}));
tally=0;
for i=1:length(alltheta_trialAv_LED.low.S)
    if ~isnan(alltheta_trialAv_LED.low.S{i})
        sumS=sumS+alltheta_trialAv_LED.low.S{i};
        tally=tally+1;
    end
end
sumS=sumS./tally;
figure(); imagesc(alltheta_trialAv_LED.t,alltheta_trialAv_LED.f(alltheta_trialAv_LED.f<=50),sumS(:,alltheta_trialAv_LED.f<=50)');
title('alltheta_trialAv_LED ALL F1');
save([outputDir '\' 'alltheta_trialAv_LED ALL F1.mat'],'sumS');
save([outputDir '\' 'UbyUalltheta_trialAv_LED.mat'],'alltheta_trialAv_LED');
% HFa_lowF1=nanmean(sumS(:,alltheta_trialAv_LED.f>=11 & alltheta_trialAv_LED.f<=20),2);
% sumS=zeros(size(alltheta_trialAv_LED.high.S{1}));
% tally=0;
% for i=1:length(alltheta_trialAv_LED.high.S)
%     if ~isnan(alltheta_trialAv_LED.high.S{i})
%         sumS=sumS+alltheta_trialAv_LED.high.S{i};
%         tally=tally+1;
%     end
% end
% sumS=sumS./tally;
% figure(); imagesc(alltheta_trialAv_LED.t,alltheta_trialAv_LED.f(alltheta_trialAv_LED.f<=50),sumS(:,alltheta_trialAv_LED.f<=30)');
% title('alltheta_trialAv_LED HIGH F1');
% save([outputDir '\' 'alltheta_trialAv_LED HIGH F1.mat'],'sumS');
% HFa_highF1=nanmean(sumS(:,alltheta_trialAv_LED.f>=11 & alltheta_trialAv_LED.f<=20),2);
% figure(); plot(alltheta_trialAv_LED.t,HFa_lowF1,'Color','k'); hold on; plot(alltheta_trialAv_LED.t,HFa_highF1,'Color','r');


end


function [lowF1_noTheta_trialAv_noLED,highF1_noTheta_trialAv_noLED]=sortLowAndHighF1(noTheta_trialAv_noLED,outputDir,tag)

F1range=[2.5 3.5];
stimWindow=[4 6.5];

% Sort into low and high F1 trials for each stimulus
for j=1:length(noTheta_trialAv_noLED)
    lowF1_noTheta_trialAv_noLED(j).t=noTheta_trialAv_noLED(j).allS.t;
    highF1_noTheta_trialAv_noLED(j).t=noTheta_trialAv_noLED(j).allS.t;
    lowF1_noTheta_trialAv_noLED(j).f=noTheta_trialAv_noLED(j).allS.f;
    highF1_noTheta_trialAv_noLED(j).f=noTheta_trialAv_noLED(j).allS.f;
    lowF1(j).t=noTheta_trialAv_noLED(j).allS.t;
    highF1(j).t=noTheta_trialAv_noLED(j).allS.t;
    lowF1(j).f=noTheta_trialAv_noLED(j).allS.f;
    highF1(j).f=noTheta_trialAv_noLED(j).allS.f;
    currt=noTheta_trialAv_noLED(1).allS.t;
    currf=noTheta_trialAv_noLED(1).allS.f;
    F1s=nan(length(noTheta_trialAv_noLED(j).allS.S),size(noTheta_trialAv_noLED(j).allS.S{1},3));
    allFs=nan(length(noTheta_trialAv_noLED(j).allS.S),size(noTheta_trialAv_noLED(j).allS.S{1},3));
    for k=1:length(noTheta_trialAv_noLED(j).allS.S)
        currs=noTheta_trialAv_noLED(j).allS.S{k};
%         if all(isnan(currs))
%             lowF1_noTheta_trialAv_noLED(j).S{k}=nan(size(currs));
%             highF1_noTheta_trialAv_noLED(j).S{k}=nan(size(currs));
%             lowF1(j).S{k}=nan(size(currs,1),size(currs,2));
%             highF1(j).S{k}=nan(size(currs,1),size(currs,2));
%             continue
%         end
        F1s(k,:)=reshape(nanmean(nanmean(currs(currt>=stimWindow(1) & currt<=stimWindow(2),currf>=F1range(1) & currf<=F1range(2),:),2),1),1,size(currs,3));
%         medF1=median(F1s);
        allFs(k,:)=reshape(nanmean(nanmean(currs(currt>=stimWindow(1) & currt<=stimWindow(2),currf>=0 & currf<=50,:),2),1),1,size(currs,3));
    end
    F1s=nanmean(F1s,1)./nanmean(allFs,1);
    medF1=median(F1s(~isnan(F1s)));
    for k=1:length(noTheta_trialAv_noLED(j).allS.S)
        currs=noTheta_trialAv_noLED(j).allS.S{k};
        if all(isnan(currs))
            lowF1_noTheta_trialAv_noLED(j).S{k}=nan(size(currs));
            highF1_noTheta_trialAv_noLED(j).S{k}=nan(size(currs));
            lowF1(j).S{k}=nan(size(currs,1),size(currs,2));
            highF1(j).S{k}=nan(size(currs,1),size(currs,2));
            continue
        end
        lowF1_noTheta_trialAv_noLED(j).S{k}=reshape(nanmean(currs(:,:,F1s<=medF1),3),size(currs,1),size(currs,2));
        highF1_noTheta_trialAv_noLED(j).S{k}=reshape(nanmean(currs(:,:,F1s>medF1),3),size(currs,1),size(currs,2));
        lowF1(j).S{k}=currs(:,:,F1s<=medF1);
        highF1(j).S{k}=currs(:,:,F1s>medF1);
    end
end
save([outputDir '\' tag 'lowF1.mat'],'lowF1');
save([outputDir '\' tag 'highF1.mat'],'highF1');


end


function avSpecgram=avS(temp_struct)

if length(size(temp_struct(1).allS.S{1}))>2
    % Need to average because didn't divide into low and high F1
    for i=1:length(temp_struct)
        for j=1:length(temp_struct(i).allS.S)
            temp_struct(i).allS.S{j}=reshape(nanmean(temp_struct(i).allS.S{j},3),size(temp_struct(i).allS.S{j},1),size(temp_struct(i).allS.S{j},2));
        end
    end
end

avSpecgram.t=[];
avSpecgram.f=[];
avSpecgram.S=[];
if length(temp_struct)==1
    avSpecgram.t=temp_struct.allS.t;
    avSpecgram.f=temp_struct.allS.f;
    avSpecgram.S=temp_struct.allS.S;
else
    avSpecgram.t=temp_struct(1).allS.t;
    avSpecgram.f=temp_struct(1).allS.f;
    avSpecgram.S=cell(1,length(temp_struct(1).allS.S));
    for i=1:length(temp_struct(1).allS.S)
        avSpecgram.S{i}=zeros(size(temp_struct(1).allS.S{1}));
    end
    for i=1:length(temp_struct)
        for j=1:length(temp_struct(1).allS.S)
            avSpecgram.S{j}=avSpecgram.S{j}+temp_struct(i).allS.S{j};
        end
    end
    for i=1:length(avSpecgram.S)
        avSpecgram.S{i}=avSpecgram.S{i}./length(temp_struct);
    end
end

end


