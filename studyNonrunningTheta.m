function studyNonrunningTheta(datadir,whichFreqInds,freqs,t,stimWindow,spontWindow,trodeDirs,runningTrials)

count_notheta=1;
count_theta=1;
count_nonrunningtheta=1;
for i=1:length(whichFreqInds)
    w=whichFreqInds(i);
%     a=load([datadir '\F1acrosscells_Hz' num2str(freqs(w))]);
    a=load([datadir '\F1acrosscells']);
    typeOfTrial='noTheta_noLED_allS_S';
    temp=a.F1acrosscells.(typeOfTrial);
    if i==1
        evF1_noTheta=temp;
    else
        evF1_noTheta=sum2arrays(evF1_noTheta,temp);
        if any(~isnan(temp(1:end)))
            count_notheta=count_notheta+1;
        end
    end
%     a=load([datadir '\F1acrosscells_Hz' num2str(freqs(w)) 'runtheta']);
%     typeOfTrial='running_theta';
    typeOfTrial='theta_noLED_allS_S';
    temp=a.F1acrosscells.(typeOfTrial);
    if i==1
        evF1_theta=temp;
    else
        evF1_theta=sum2arrays(evF1_theta,temp);
        if any(~isnan(temp(1:end)))
            count_theta=count_theta+1;
        end
    end
%     a=load([datadir '\F1acrosscells_Hz' num2str(freqs(w))]);
    typeOfTrial='nonrunning_theta';
    temp=a.F1acrosscells.(typeOfTrial);
    if i==1
        evF1_nonrunningtheta=temp;
    else
        evF1_nonrunningtheta=sum2arrays(evF1_nonrunningtheta,temp);
        if any(~isnan(temp(1:end)))
            count_nonrunningtheta=count_nonrunningtheta+1;
        end
    end
end
evF1_noTheta=evF1_noTheta/count_notheta;
evF1_theta=evF1_theta/count_theta;
evF1_nonrunningtheta=evF1_nonrunningtheta/count_nonrunningtheta;

temp1=nanmean(evF1_noTheta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_noTheta(:,t>spontWindow(1) & t<spontWindow(2)),2);
temp2=nanmean(evF1_theta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_theta(:,t>spontWindow(1) & t<spontWindow(2)),2);
temp3=nanmean(evF1_nonrunningtheta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_nonrunningtheta(:,t>spontWindow(1) & t<spontWindow(2)),2);
% takeUnits=temp1>0 | temp2>0 | temp3>0;
% takeUnits=temp1>0 & temp2>0 & temp3>0;
takeUnits=1:size(evF1_noTheta,1);
figure(); plot(t,nanmean(evF1_noTheta(takeUnits,:),1),'Color','k'); hold on;
plot(t,nanmean(evF1_noTheta(takeUnits,:),1)-nanstd(evF1_noTheta(takeUnits,:),[],1)./sqrt(size(evF1_noTheta(takeUnits,:),1)),'Color','k');
plot(t,nanmean(evF1_noTheta(takeUnits,:),1)+nanstd(evF1_noTheta(takeUnits,:),[],1)./sqrt(size(evF1_noTheta(takeUnits,:),1)),'Color','k');
plot(t,nanmean(evF1_theta(takeUnits,:),1),'Color','r');
plot(t,nanmean(evF1_theta(takeUnits,:),1)-nanstd(evF1_theta(takeUnits,:),[],1)./sqrt(size(evF1_theta(takeUnits,:),1)),'Color','r');
plot(t,nanmean(evF1_theta(takeUnits,:),1)+nanstd(evF1_theta(takeUnits,:),[],1)./sqrt(size(evF1_theta(takeUnits,:),1)),'Color','r');
plot(t,nanmean(evF1_nonrunningtheta(takeUnits,:),1),'Color','b');
plot(t,nanmean(evF1_nonrunningtheta(takeUnits,:),1)-nanstd(evF1_nonrunningtheta(takeUnits,:),[],1)./sqrt(size(evF1_nonrunningtheta(takeUnits,:),1)),'Color','b');
plot(t,nanmean(evF1_nonrunningtheta(takeUnits,:),1)+nanstd(evF1_nonrunningtheta(takeUnits,:),[],1)./sqrt(size(evF1_nonrunningtheta(takeUnits,:),1)),'Color','b');

temp1=nanmean(evF1_noTheta(takeUnits,:),1);
temp2=nanmean(evF1_theta(takeUnits,:),1);
temp3=nanmean(evF1_nonrunningtheta(takeUnits,:),1);
figure(); 
temp1=nanmean(temp1(t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(temp1(t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp1;
scatter(1,nanmean(temp),[],'k');
hold on;
temp2=nanmean(temp2(t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(temp2(t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp2;
scatter(3,nanmean(temp),[],'r');
temp3=nanmean(temp3(t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(temp3(t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp3;
scatter(2,nanmean(temp),[],'b');
line([1 2 3],[nanmean(temp1) nanmean(temp3) nanmean(temp2)]);
title('mean across cells first evoked');

figure(); 
temp1=nanmean(evF1_noTheta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_noTheta(:,t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp1;
scatter(1,nanmean(temp),[],'k');
line([1 1],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','k');
hold on;
temp2=nanmean(evF1_theta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_theta(:,t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp2;
scatter(3,nanmean(temp),[],'r');
line([3 3],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','r');
temp3=nanmean(evF1_nonrunningtheta(:,t>stimWindow(1) & t<stimWindow(2)),2)-nanmean(evF1_nonrunningtheta(:,t>spontWindow(1) & t<spontWindow(2)),2);
temp=temp3;
scatter(2,nanmean(temp),[],'b');
line([2 2],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','b');
line([1 2 3],[nanmean(temp1) nanmean(temp3) nanmean(temp2)]);
title('evoked');

noTheta_spont=[];
theta_spont=[];
nonrunningtheta_spont=[];
noTheta_psths=[];
theta_psths=[];
nonrunningtheta_psths=[];
for i=1:length(trodeDirs)
    a=load([datadir '\' trodeDirs{i} '\dLGNpsth.mat']);
    dLGNpsth=a.dLGNpsth;
    a=load([datadir '\' trodeDirs{i} '\noThetaTrials.mat']);
    noThetaTrials=a.noThetaTrials;
    for j=1:length(dLGNpsth.psths)
        l=dLGNpsth.unitLED{1};
        temp=dLGNpsth.psths{j};
        noTheta_psths=[noTheta_psths; nanmean(temp(ismember(l,freqs) & noThetaTrials'==1,:),1)];
        theta_psths=[theta_psths; nanmean(temp(ismember(l,freqs)  & noThetaTrials'==0 & runningTrials==1,:),1)];
        nonrunningtheta_psths=[nonrunningtheta_psths; nanmean(temp(ismember(l,freqs)  & noThetaTrials'==0 & runningTrials==0,:),1)];
        noTheta_spont=[noTheta_spont; nanmean(nanmean(temp(ismember(l,freqs) & noThetaTrials'==1,:),1))];
        theta_spont=[theta_spont; nanmean(nanmean(temp(ismember(l,freqs)  & noThetaTrials'==0 & runningTrials==1,:),1))];
        nonrunningtheta_spont=[nonrunningtheta_spont; nanmean(nanmean(temp(ismember(l,freqs)  & noThetaTrials'==0 & runningTrials==0,:),1))];
    end
end
figure();
plot(dLGNpsth.t,nanmean(noTheta_psths,1),'Color','k');
hold on;
plot(dLGNpsth.t,nanmean(noTheta_psths,1)-nanstd(noTheta_psths,[],1)./sqrt(size(noTheta_psths,1)),'Color','k');
plot(dLGNpsth.t,nanmean(noTheta_psths,1)+nanstd(noTheta_psths,[],1)./sqrt(size(noTheta_psths,1)),'Color','k');

plot(dLGNpsth.t,nanmean(theta_psths,1),'Color','r');
plot(dLGNpsth.t,nanmean(theta_psths,1)-nanstd(theta_psths,[],1)./sqrt(size(theta_psths,1)),'Color','r');
plot(dLGNpsth.t,nanmean(theta_psths,1)+nanstd(theta_psths,[],1)./sqrt(size(theta_psths,1)),'Color','r');

plot(dLGNpsth.t,nanmean(nonrunningtheta_psths,1),'Color','b');
plot(dLGNpsth.t,nanmean(nonrunningtheta_psths,1)-nanstd(nonrunningtheta_psths,[],1)./sqrt(size(nonrunningtheta_psths,1)),'Color','b');
plot(dLGNpsth.t,nanmean(nonrunningtheta_psths,1)+nanstd(nonrunningtheta_psths,[],1)./sqrt(size(nonrunningtheta_psths,1)),'Color','b');

figure();
temp1=noTheta_spont;
temp=temp1;
scatter(1,nanmean(temp),[],'k');
line([1 1],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','k');
hold on;
temp2=theta_spont;
temp=temp2;
scatter(3,nanmean(temp),[],'r');
line([3 3],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','r');
temp3=nonrunningtheta_spont;
temp=temp3;
scatter(2,nanmean(temp),[],'b');
line([2 2],[nanmean(temp)-nanstd(temp)./length(temp) nanmean(temp)+nanstd(temp)./length(temp)],'Color','b');
line([1 2 3],[nanmean(temp1) nanmean(temp3) nanmean(temp2)]);
title('spont');

end

function C=sum2arrays(A,B)

tmp = cat(3,A,B); 
C = nansum(tmp,3);

end
