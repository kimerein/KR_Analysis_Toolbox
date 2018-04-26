function plotSummaryOfSUFreqResponse(psth,ind,uses,uses_tri,noThetaTrials,isNoTheta)

% % ds=1;
% ds=5; % used for temporal frequency 3 Hz in paper
% ds=1;
% % smoo=1;
% smoo=3; % used for temporal frequency 3 Hz in paper
% smoo=1;
% leduse=[20.05];
% noleduse=[20.00];

figure();

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
ledfreqs=freqs+0.05;
% freqs=[0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.20 0.30 0.40 0.50 0.60];
% ledfreqs=freqs+5;

ha = tight_subplot(5,3,[.03 .05],[.05 .02],[.03 .03]);
for i=1:length(freqs)
    axes(ha(i)); 
    if ismember(i,[1 2])
        plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,5,3,ledfreqs(i),freqs(i),1.5);
    elseif ismember(i,[3 4])
%         plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,3,2,ledfreqs(i),freqs(i),1.5);
        plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,5,2,ledfreqs(i),freqs(i),1.5);
    elseif ismember(i,[5 6])
        plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,2,1,ledfreqs(i),freqs(i),1.5);
    else
%         plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,1,1,ledfreqs(i),freqs(i),1.5);
        plotSUpsth(psth,uses,uses_tri,ind,noThetaTrials,isNoTheta,1,2,ledfreqs(i),freqs(i),1.5);
    end
end
set(ha([1:12]),'XTickLabel',''); 
% set(ha,'YTickLabel','');

end

function [black_t,red_t,black_av,red_av]=plotSUpsth(dLGNpsth,uses,uses_tri,i,noThetaTrials,isNoTheta,ds,smoo,leduse,noleduse,lw)

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
%     figure();
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(l==noleduse & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','k');
%     hold on;
%     plot(smooth(downSampAv(t,ds),smoo),smooth(downSampAv(nanmean(p(ismember(l,leduse) & noThetaTrials'==isNoTheta & takeTrials,:),1),ds),smoo),'Color','b');
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

% figure();
plot(smooth(downSampAv(t,ds),smoo),runningav,'Color','k','LineWidth',lw);
black_t=smooth(downSampAv(t,ds),smoo);
black_av=runningav;
hold on;
plot(smooth(downSampAv(t,ds),smoo),runningav_led,'Color',[0 0.5 1],'LineWidth',lw);
red_t=smooth(downSampAv(t,ds),smoo);
red_av=runningav_led;

% plot(smooth(downSampAv(t,ds),smoo),runningav,'Color','k','LineWidth',lw);
% black_t=smooth(downSampAv(t,ds),smoo);
% black_av=runningav;
xlim([0 4]);

end