function KER_plotOrganizedLFP(LFPdata,Fs,redTrials,blueTrials)

trialLength=size(LFPdata,2)*(1/Fs);

% Align all LFP traces to initial value = 0
for i=1:size(LFPdata,1)
    initVal=LFPdata(i,1);
    LFPdata(i,:)=LFPdata(i,:)-initVal;
end

figure;
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata,1),'Color','k');
title('Average Stimulus-Triggered LFP - Initial Values Aligned');

% Plot averaged during diff. brain states with std. devs
figure;
redTrials=redTrials';
blueTrials=blueTrials';
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:),1),'Color','r');
hold on;
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:),1)+std(LFPdata(redTrials,:),1),'Color','m');
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:),1)-std(LFPdata(redTrials,:),1),'Color','m');
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:),1)+std(LFPdata(blueTrials,:),1),'Color','c');
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:),1)-std(LFPdata(blueTrials,:),1),'Color','c');

figure;
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(redTrials,:),1),'Color','r');
hold on;
plot(0:trialLength/(size(LFPdata,2)-1):trialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
minValForPlot=min(min([mean(LFPdata(redTrials,:),1); mean(LFPdata(blueTrials,:),1)],[],2),[],1);
for i=0.02:0.02:trialLength
    [h,p,ci]=ttest2(mean(LFPdata(redTrials,floor((i-0.02)/(trialLength/(size(LFPdata,2)-1)))+1:floor(i/(trialLength/(size(LFPdata,2)-1)))+1),2),mean(LFPdata(blueTrials,floor((i-0.02)/(trialLength/(size(LFPdata,2)-1)))+1:floor(i/(trialLength/(size(LFPdata,2)-1)))+1),2));
    if h==1
        line([i-0.02 i],[minValForPlot minValForPlot],'Color','k');
    end
end
