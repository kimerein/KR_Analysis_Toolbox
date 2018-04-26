function [black_t,red_t,black_av,red_av]=plotSUpsth(dLGNpsth,uses,uses_tri,i,noThetaTrials,isNoTheta)

% ds=1;
ds=5; % used for temporal frequency 3 Hz in paper
ds=1;
% smoo=1;
smoo=3; % used for temporal frequency 3 Hz in paper
smoo=1;
leduse=[20.05];
noleduse=[20.00];

t=dLGNpsth.t;
l=dLGNpsth.unitLED{1};
p=dLGNpsth.psths{i};
tri=dLGNpsth.unitTrials{1};
s=dLGNpsth.unitStimcond{1};

runningav=[];
runningav_led=[];
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    takeTrials=zeros(1,length(tri));
    for j=1:length(currs)
        currcurrs=currs{j};
        currcurrtri=currtri{j};
        takeTrials(ismember(s,currcurrs) & ismember(tri,currcurrtri))=1;
    end
    figure();
    plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==noleduse & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','k');
    hold on;
    plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(ismember(l,leduse) & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','b');
    if i==1
        runningav_led=smooth(downSampAv(nanmean(p(ismember(l,leduse) & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
        runningav=smooth(downSampAv(nanmean(p(l==noleduse & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
    else
        runningav_led=runningav_led+smooth(downSampAv(nanmean(p(ismember(l,leduse) & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
        runningav=runningav+smooth(downSampAv(nanmean(p(l==noleduse & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
    end
end
runningav=runningav./length(uses);
runningav_led=runningav_led./length(uses);

figure();
plot(smooth(downSampAv(t,ds),smoo),runningav,'Color','k');
black_t=smooth(downSampAv(t,ds),smoo);
black_av=runningav;
hold on;
plot(smooth(downSampAv(t,ds),smoo),runningav_led,'Color','b');
red_t=smooth(downSampAv(t,ds),smoo);
red_av=runningav_led;
title('Average');