function [black_t,red_t,black_av,red_av]=plotSUpsth_withPref(dLGNpsth,uses,uses_tri,i,noThetaTrials,isNoTheta,allspikespsth)

ds=5;
% ds=10;
smoo=3;

t=dLGNpsth.t;
l=dLGNpsth.unitLED{1};
p=dLGNpsth.psths{i};
p_for_pref=allspikespsth.psths{i};
tri=dLGNpsth.unitTrials{1};
s=dLGNpsth.unitStimcond{1};

runningav=[];
runningav_led=[];
for i=1:length(uses)
    currs=uses{i};
    currtri=uses_tri{i};
    takeTrials=zeros(1,length(tri));
    psths_for_pref=[];
    for j=1:length(currs)
        currcurrs=currs{j};
        currcurrtri=currtri{j};
        takeTrials(ismember(s,currcurrs) & ismember(tri,currcurrtri))=1;
    end
    psths_for_pref(i,:)=nanmean(p_for_pref(l==0 & noThetaTrials'==0 & takeTrials,:),1);
%     figure();
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==0 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','k');
%     hold on;
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==5.05 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','b');
    if i==1
        runningav_led=smooth(downSampAv(nanmean(p(l==5.05 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
        runningav=smooth(downSampAv(nanmean(p(l==0 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
    else
        runningav_led=runningav_led+smooth(downSampAv(nanmean(p(l==5.05 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
        runningav=runningav+smooth(downSampAv(nanmean(p(l==0 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
    end
end

responses=nanmean(psths_for_pref(:,t>=4 & t<=6.5),2)-nanmean(psths_for_pref(:,t>=0 & t<=4),2);
[ma,mind]=max(responses);


runningav=[];
runningav_led=[];
for i=mind
    currs=uses{i};
    currtri=uses_tri{i};
    takeTrials=zeros(1,length(tri));
%     psths_for_pref=[];
    for j=1:length(currs)
        currcurrs=currs{j};
        currcurrtri=currtri{j};
        takeTrials(ismember(s,currcurrs) & ismember(tri,currcurrtri))=1;
    end
%     psths_for_pref(i,:)=nanmean(p_for_pref(l==0 & noThetaTrials'==0 & takeTrials,:),1)
%     figure();
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==0 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','k');
%     hold on;
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==5.05 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','b');
    runningav_led=smooth(downSampAv(nanmean(p(l==5.05 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
    runningav=smooth(downSampAv(nanmean(p(l==0 & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo);
end


figure();
plot(smooth(downSampAv(t,ds),smoo),runningav,'Color','k');
black_t=smooth(downSampAv(t,ds),smoo);
black_av=runningav;
hold on;
plot(smooth(downSampAv(t,ds),smoo),runningav_led,'Color','b');
red_t=smooth(downSampAv(t,ds),smoo);
red_av=runningav_led;