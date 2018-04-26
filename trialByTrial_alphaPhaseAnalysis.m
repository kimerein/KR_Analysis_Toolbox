function [currTrials,nextTrials,prevTrials,bootout,outbins,outps,outps_prev,outsum]=trialByTrial_alphaPhaseAnalysis(prefStimset,psths_trialByTrial,allRef_trialByTrial,sigStim,arrayVec)

xt=psths_trialByTrial.xt;
psths_t=psths_trialByTrial.psths_t;
psths_ref=allRef_trialByTrial.psths_t;

% Show example trials of spiking in which you have separated units' spikes by
% tuning to stimuli
% TUNING to PHASE

% Re-organize psths by current stim.
curr_psths=cell(size(psths_t,1),1);
x1=xt(1):xt(2)-xt(1):4.5;
y=nan(size(psths_t{1,1},1)+size(psths_t{1,2},1)+size(psths_t{1,3},1),length(xt(xt>=0 & xt<=1.5))*3);
allRef_psths=nan(size(psths_t{1,1},1)+size(psths_t{1,2},1)+size(psths_t{1,3},1),length(xt(xt>=0 & xt<=1.5))*3);
for i=1:size(psths_t,1)
    k=1;
    curry=psths_t{i,1};
    if i==1 && ~isempty(allRef_trialByTrial.xt)
        curry_ref=psths_ref{1};
        tempref=makeArraySameLengthAsVec(x1(x1>=0 & x1<=3.5),curry_ref);
        allRef_psths(k:k+size(tempref,1)-1,x1>=0 & x1<=3.5)=tempref;
    end
    temp=makeArraySameLengthAsVec(x1(x1>=0 & x1<=3.5),curry);
    y(k:k+size(temp,1)-1,x1>=0 & x1<=3.5)=temp;
    k=k+size(temp,1);
    curry=psths_t{i,2};
    if i==1 && ~isempty(allRef_trialByTrial.xt)
        curry_ref=psths_ref{2};
        tempref=makeArraySameLengthAsVec(x1(x1>=0 & x1<=0.5),curry_ref(:,xt>=3));
        allRef_psths(k:k+size(tempref,1)-1,x1>=0 & x1<=0.5)=tempref;
        tempref=makeArraySameLengthAsVec(x1(x1>=1.5 & x1<=4.5),curry_ref(:,xt>=0 & xt<=3));
        allRef_psths(k:k+size(tempref,1)-1,x1>=1.5 & x1<=4.5)=tempref;
    end
    temp=makeArraySameLengthAsVec(x1(x1>=0 & x1<=0.5),curry(:,xt>=3));
    y(k:k+size(temp,1)-1,x1>=0 & x1<=0.5)=temp;
    curry=psths_t{i,2};
    temp=makeArraySameLengthAsVec(x1(x1>=1.5 & x1<=4.5),curry(:,xt>=0 & xt<=3));
    y(k:k+size(temp,1)-1,x1>=1.5 & x1<=4.5)=temp;
    k=k+size(temp,1);   
    if i==1 && ~isempty(allRef_trialByTrial.xt)
        curry_ref=psths_ref{3};
        tempref=makeArraySameLengthAsVec(x1(x1>=0 & x1<=2),curry_ref(:,xt>=1.5 & xt<=3.5));
        allRef_psths(k:k+size(tempref,1)-1,x1>=0 & x1<=2)=tempref;
        curry_ref=psths_ref{3};
        tempref=makeArraySameLengthAsVec(x1(x1>=3 & x1<=4.5),curry_ref(:,xt>=0 & xt<=1.5));
        allRef_psths(k:k+size(tempref,1)-1,x1>=3 & x1<=4.5)=tempref;
    end
    curry=psths_t{i,3};
    temp=makeArraySameLengthAsVec(x1(x1>=0 & x1<=2),curry(:,xt>=1.5 & xt<=3.5));
    y(k:k+size(temp,1)-1,x1>=0 & x1<=2)=temp;
    curry=psths_t{i,3};
    temp=makeArraySameLengthAsVec(x1(x1>=3 & x1<=4.5),curry(:,xt>=0 & xt<=1.5));
    y(k:k+size(temp,1)-1,x1>=3 & x1<=4.5)=temp;
    curr_psths{i}=y;
end

% getPvaluePSTHforEachUnit(curr_psths,prefStimset,sigStim,x1,arrayVec);
getPvalueTrials(curr_psths,prefStimset,sigStim,x1,arrayVec);

% Add up spikes from units with same tuning
% tuned1trials=addUnitsByTuning(curr_psths,prefStimset,1,sigStim,arrayVec);
% tuned2trials=addUnitsByTuning(curr_psths,prefStimset,2,sigStim,arrayVec);
% tuned3trials=addUnitsByTuning(curr_psths,prefStimset,3,sigStim,arrayVec);
tuned1trials=addUnitsByTuning(curr_psths,prefStimset,1,sigStim,[]);
tuned2trials=addUnitsByTuning(curr_psths,prefStimset,2,sigStim,[]);
tuned3trials=addUnitsByTuning(curr_psths,prefStimset,3,sigStim,[]);

% Bootstrap weighted representation
bootout=bootstrapWeighted(curr_psths,arrayVec,x1);

% Bootstrap sum
outsum=bootstrapSum(curr_psths,prefStimset,sigStim,x1,arrayVec);

% Reorganize data according to current, previous and next
currTrials=tuned1trials(:,x1>=0 & x1<=1.5);
currTrials=[currTrials; tuned2trials(:,x1>=1.5 & x1<=3)];
currTrials=[currTrials; tuned3trials(:,x1>=3 & x1<=4.5)];
nextTrials=tuned3trials(:,x1>=0 & x1<=1.5);
nextTrials=[nextTrials; tuned1trials(:,x1>=1.5 & x1<=3)];
nextTrials=[nextTrials; tuned2trials(:,x1>=3 & x1<=4.5)];
prevTrials=tuned2trials(:,x1>=0 & x1<=1.5);
prevTrials=[prevTrials; tuned3trials(:,x1>=1.5 & x1<=3)];
prevTrials=[prevTrials; tuned1trials(:,x1>=3 & x1<=4.5)];
refTrials=allRef_psths(:,x1>=0 & x1<=1.5);
refTrials=[refTrials; allRef_psths(:,x1>=1.5 & x1<=3)];
refTrials=[refTrials; allRef_psths(:,x1>=3 & x1<=4.5)];

% nextTrials=prevTrials;
% Plot some example trials
d=1;
newCurrTrials=zeros(size(currTrials,1),length(downSampAv(currTrials(1,:),d)));
newNextTrials=zeros(size(currTrials,1),length(downSampAv(currTrials(1,:),d)));
newx=downSampAv(x1,d);
newCurrTrials=downSampMatrix(currTrials,d);
newNextTrials=downSampMatrix(nextTrials,d);
newPrevTrials=downSampMatrix(prevTrials,d);
figure();
plot(newx(newx>=1.0 & newx<=1.5),newNextTrials(:,newx>=1.0 & newx<=1.5)'+rand(size(newNextTrials(:,newx>=1.0 & newx<=1.5)))'.*20,'Color','r');
hold on;
plot(newx(newx>=1.0 & newx<=1.5),newCurrTrials(:,newx>=1.0 & newx<=1.5)'+rand(size(newCurrTrials(:,newx>=1.0 & newx<=1.5)))'.*20,'Color','k');
% scatter(newx(newx>=1.2 & newx<=1.5),newNextTrials(:,newx>=1.2 & newx<=1.5)'+rand(size(newNextTrials(:,newx>=1.2 & newx<=1.5)))'.*20,[],'r');
% hold on;
% scatter(newx(newx>=1.2 & newx<=1.5),newCurrTrials(:,newx>=1.2 & newx<=1.5)'+rand(size(newCurrTrials(:,newx>=1.2 & newx<=1.5)))'.*20,[],'k');
% hold on;

figure(); 
plot(newx(newx>=1.0 & newx<=1.5),nansum(newNextTrials(:,newx>=1.0 & newx<=1.5),1),'Color','r');
hold on;
plot(newx(newx>=1.0 & newx<=1.5),nansum(newCurrTrials(:,newx>=1.0 & newx<=1.5),1),'Color','k');

offsets=[0 0.05 0.1 0.15];
tryBins=sort(1.5:-0.2:0);
allp=nan(length(offsets),length(tryBins)-1);
o=cell(1,length(offsets));
allp_prev=nan(length(offsets),length(tryBins)-1);
for j=1:length(offsets)
    tryBins=sort(1.5-offsets(j):-0.2:0);
    o{j}=tryBins;
    % tryBins=0:0.167:1.5;
    pvalsAcrossTrial=ones(1,length(tryBins)-1);
    pvalsAcrossTrial_prev=ones(1,length(tryBins)-1);
    for i=1:length(tryBins)-1
        window=[tryBins(i) tryBins(i+1)];
        pvalsAcrossTrial(i)=compareDistributions(newx,newNextTrials,newCurrTrials,window);
        pvalsAcrossTrial_prev(i)=compareDistributions(newx,newPrevTrials,newCurrTrials,window);
    end
    allp(j,end-length(pvalsAcrossTrial)+1:end)=pvalsAcrossTrial;
    allp_prev(j,end-length(pvalsAcrossTrial_prev)+1:end)=pvalsAcrossTrial_prev;
end
figure(); 
hax=axes();
outbins=[];
outps=[];
outps_prev=[];
for i=1:size(allp,1)
    oo=o{i};
    scatter(nanmean([oo(1:end-1); oo(2:end)],1),allp(i,end-length(oo(1:end-1))+1:end)+rand(size(allp(i,end-length(oo(1:end-1))+1:end))).*0.01,[],'k');
    outbins=[outbins nanmean([oo(1:end-1); oo(2:end)],1)];
    outps=[outps allp(i,end-length(oo(1:end-1))+1:end)+rand(size(allp(i,end-length(oo(1:end-1))+1:end))).*0.01];
    outps_prev=[outps_prev allp_prev(i,end-length(oo(1:end-1))+1:end)+rand(size(allp_prev(i,end-length(oo(1:end-1))+1:end))).*0.01];
    hold on;
end
% plot(mean([tryBins(1:end-1); tryBins(2:end)],1),pvalsAcrossTrial);
title('pvals across trial');
disp(allp);

nexamp=1;
useTrials1=1+floor(rand(1,nexamp).*(size(currTrials,1)-1));
k=1;
while ~any(currTrials(useTrials1,x1>=1.2 & x1<=1.5)~=0 & ~isnan(currTrials(useTrials1,x1>=1.2 & x1<=1.5)))
    useTrials1=1+floor(rand(1,nexamp).*(size(currTrials,1)-1));
    k=k+1;
    if k>10000
        break
    end
end
useTrials=useTrials1;
figure(); 
plot(downSampAv(x1(x1>=1.2 & x1<=1.5),5),downSampAv(currTrials(useTrials,x1>=1.2 & x1<=1.5),5)','Color','k');
hold on;
plot(downSampAv(x1(x1>=1.2 & x1<=1.5),5),downSampAv(nextTrials(useTrials,x1>=1.2 & x1<=1.5),5)','Color','r');
% plot(x1(x1>=0 & x1<=1.5),prevTrials(useTrials,:)','Color','g');
% plot(x1(x1>=0 & x1<=1.5),nanmean(refTrials(useTrials,:),1)','Color','c');

% % Cycle-average phase, but do trial-by-trial, then align across trials
% params.Fs=1/(x1(2)-x1(1));
% params.tapers=[0.5 x1(end)-x1(1) 0];
% [avSp,avt,avf]=mtspecgrampb(nanmean(currTrials,1),1),[0.5 0.1],params); % what is this 0.5 0.1 thing?
% allSp=nan(size(currTrials,1),length(avSp));
% allf=nan(size(currTrials,1),length(avf));
% allt=nan(size(currTrials,1),length(avt));
% for i=1:size(currTrials,1)
%       [Sp,t,f]=mtspecgrampb(currTrials(i,:));
% end
% 

end

function p=compareDistributions(newx,newNextTrials,newCurrTrials,window)

b=nansum(newNextTrials(:,newx>window(1) & newx<=window(2)),1);
a=nansum(newCurrTrials(:,newx>=window(1) & newx<=window(2)),1);
useX=newx(newx>=window(1) & newx<=window(2));
nextTrialsSet=[];
for i=1:length(b)
    addN=ceil(b(i)/60);
    for j=1:addN
        nextTrialsSet=[nextTrialsSet useX(i)];
    end
end
currTrialsSet=[];
for i=1:length(a)
    addN=ceil(a(i)/60);
    for j=1:addN
        currTrialsSet=[currTrialsSet useX(i)];
    end
end
p=ranksum(currTrialsSet,nextTrialsSet);
% disp('pvalue for curr and next distributions');
% disp(p);

end

function out=getPvalueTrials(data,prefStimset,sigStim,x1,arrayVec)

binsize=1;
out=[];
colorSigCutoff=0.05;
colorSigCutoff=1;
colorByPref=0;
chooseTrial=298;
% colorSigCutoff=0.0001;

% psth1vs2=nan(length(data),size(d,1),size(downSampMatrix(data{1},binsize),2));
% psth1vs3=nan(length(data),size(d,1),size(downSampMatrix(data{1},binsize),2));
% psth2vs3=nan(length(data),size(d,1),size(downSampMatrix(data{1},binsize),2));
psth1vsall=nan(length(data),size(data{1},1),size(downSampMatrix(data{1},binsize),2));
psth2vsall=nan(length(data),size(data{1},1),size(downSampMatrix(data{1},binsize),2));
psth3vsall=nan(length(data),size(data{1},1),size(downSampMatrix(data{1},binsize),2));
subx=downSampMatrix(x1(x1>=0 & x1<=1.5),binsize);
allData=zeros(size(data{1}));
for i=1:length(data)
    d=data{i};
    allData=allData+d;
end
fTrials=find(nansum(allData,2)>0);
fTrials=sort(fTrials);
fTrials=unique(fTrials);

for i=1:length(data)
    d=data{i};
    sub1=downSampMatrix(d(:,x1>=0 & x1<=1.5),binsize);
    sub2=downSampMatrix(d(:,x1>=1.5 & x1<=3),binsize);
    sub3=downSampMatrix(d(:,x1>=3 & x1<=4.5),binsize);
    for j=1:length(subx)
        a=sub1(:,j);
        b=sub2(:,j);
        c=sub3(:,j);
        for k=fTrials(chooseTrial)
            if isnan(a(k)) || (any(b(~isnan(b))~=0)==0 && any(c(~isnan(c))~=0)==0)
                p=1;
            else
                [~,p]=ttest2(a(k),[a(~isnan(a)); b(~isnan(b)); c(~isnan(c))]);
            end
            if a(k)>nanmean([b(~isnan(b)); c(~isnan(c))])
                psth1vsall(i,k,j)=p;
            else
                psth1vsall(i,k,j)=-p;
            end
            if isnan(b(k)) || (any(a(~isnan(a))~=0)==0 && any(c(~isnan(c))~=0)==0)
                p=1;
            else
                [~,p]=ttest2(b(k),[a(~isnan(a)); b(~isnan(b)); c(~isnan(c))]);
            end
            if b(k)>nanmean([a(~isnan(a)); c(~isnan(c))])
                psth2vsall(i,k,j)=p;
            else
                psth2vsall(i,k,j)=-p;
            end
            if isnan(c(k)) || (any(b(~isnan(b))~=0)==0 && any(a(~isnan(a))~=0)==0)
                p=1;
            else
                [~,p]=ttest2(c(k),[b(~isnan(b)); a(~isnan(a)); c(~isnan(c))]);
            end
            if c(k)>nanmean([b(~isnan(b)); a(~isnan(a))])
                psth3vsall(i,k,j)=p;
            else
                psth3vsall(i,k,j)=-p;
            end
        end
    end
end


temp1=reshape(psth1vsall(:,fTrials(chooseTrial),:),size(psth1vsall,1),size(psth1vsall,3));
temp2=reshape(psth2vsall(:,fTrials(chooseTrial),:),size(psth1vsall,1),size(psth1vsall,3));
temp3=reshape(psth3vsall(:,fTrials(chooseTrial),:),size(psth1vsall,1),size(psth1vsall,3));
clear psth1vsall psth2vsall psth3vsall
psth1vsall=temp1;
psth2vsall=temp2;
psth3vsall=temp3;

cmap1=[0 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1];
% cmap2=[1      0.9   0.8  0.7  0.6 0.5 0.4 0.3 0.2 0.1];
% cmap2=linspace(1,0.2,10);
cmap2=[1      0.95   0.9  0.85  0.8 0.75 0.7 0.65 0.6 0.55];
% cmap2=[1      1      1      1      1  1    1   1    1   1];
% Make current vs next figure
h1=figure(1); 
hax=axes();
set(gcf,'Visible','off');
title('Current (blue) vs Next (red)');
% h2=figure(2);
% hax2=axes();
% set(gcf,'Visible','off');
yoffset=0;
xoffset=subx(1)-(subx(2)-subx(1));
subx=[xoffset subx];
forimscR=zeros(length(data),length(subx)-1);
forimscB=zeros(length(data),length(subx)-1);
for i=1:length(data)
    xoffset=subx(1)-(subx(2)-subx(1));
    d=data{i};
%     sub1=nanmean(downSampMatrix(d(:,x1>=0 & x1<=1.5),binsize),1);
%     sub2=nanmean(downSampMatrix(d(:,x1>=1.5 & x1<=3),binsize),1);
%     sub3=nanmean(downSampMatrix(d(:,x1>=3 & x1<=4.5),binsize),1);
    sub1=nanmean(downSampMatrix(d(fTrials(chooseTrial),x1>=0 & x1<=1.5),binsize),1);
    sub2=nanmean(downSampMatrix(d(fTrials(chooseTrial),x1>=1.5 & x1<=3),binsize),1);
    sub3=nanmean(downSampMatrix(d(fTrials(chooseTrial),x1>=3 & x1<=4.5),binsize),1);
    maxy=max([sub1 sub2 sub3],[],2);
    % Stim 1 vs 2
    if arrayVec(i,1)>=0 & 1-arrayVec(i,1)<=colorSigCutoff % Fires more for stim 1 than 2, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,1)<=0 & 1+arrayVec(i,1)<=colorSigCutoff % Fires more for stim 2 than 1, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth1vsall(i,j)=psth1vs2(i,j);
        if psth1vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 1
            cmask=[0 0 1];
        elseif psth1vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 1 
            cmask=[0 0 1];
        elseif psth1vsall(i,j)>=0 & currSig==1 & upward==0 % If FR increase and increase means stim 2
            cmask=[1 0 0];
        elseif psth1vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 2
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        if isnan(psth1vsall(i,j))
            psth1vsall(i,j)=1;
        end
        colorsatp=abs(psth1vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if colorByPref==1
            if prefStimset(i)==1
                currc=[1 0 0];
            elseif prefStimset(i)==2
                currc=[0 1 0];
            elseif prefStimset(i)==3
                currc=[0 0 1];
            end
        end
        if j==1
            figure(1); 
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub1(j) sub1(j)],'Color',currc,'LineWidth',1.5);
            hold on;
%             figure(2);
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
%             hold on;
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else 
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub1(j-1) sub1(j)],'Color',currc,'LineWidth',1.5); 
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    % Stim 2 vs 3
    xoffset=1.5;
    if arrayVec(i,3)>=0 & 1-arrayVec(i,3)<=colorSigCutoff % Fires more for stim 2 than 3, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,3)<=0 & 1+arrayVec(i,3)<=colorSigCutoff % Fires more for stim 3 than 2, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth2vsall(i,j)=psth2vs3(i,j);
        if psth2vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 2
            cmask=[0 0 1];
        elseif psth2vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 2 
            cmask=[0 0 1];
        elseif psth2vsall(i,j)>=0 & currSig==1 & upward==0 % If FR increase and increase means stim 3
            cmask=[1 0 0];
        elseif psth2vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 3
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        if isnan(psth2vsall(i,j))
            psth2vsall(i,j)=1;
        end
        colorsatp=abs(psth2vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if colorByPref==1
            if prefStimset(i)==1
                currc=[1 0 0];
            elseif prefStimset(i)==2
                currc=[0 1 0];
            elseif prefStimset(i)==3
                currc=[0 0 1];
            end
        end
        if j==1
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub2(j) sub2(j)],'Color',currc,'LineWidth',1.5);
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub2(j-1) sub2(j)],'Color',currc,'LineWidth',1.5); 
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    % Stim 1 vs 3 (or 3 vs 1, correct order)
    xoffset=3;
    if arrayVec(i,2)>=0 & 1-arrayVec(i,2)<=colorSigCutoff % Fires more for stim 1 than 3, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,2)<=0 & 1+arrayVec(i,2)<=colorSigCutoff % Fires more for stim 3 than 1, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth3vsall(i,j)=-psth1vs3(i,j);
        if psth3vsall(i,j)>=0 & currSig==1 & upward==0 % If FR relative increase during 3 and increase means stim 3
            cmask=[0 0 1];
        elseif psth3vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 3 
            cmask=[0 0 1];
        elseif psth3vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 1
            cmask=[1 0 0];
        elseif psth3vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 1
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        if isnan(psth3vsall(i,j))
            psth3vsall(i,j)=1;
        end
        colorsatp=abs(psth3vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if colorByPref==1
            if prefStimset(i)==1
                currc=[1 0 0];
            elseif prefStimset(i)==2
                currc=[0 1 0];
            elseif prefStimset(i)==3
                currc=[0 0 1];
            end
        end
        if j==1
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub3(j) sub3(j)],'Color',currc,'LineWidth',1.5);
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub3(j-1) sub3(j)],'Color',currc,'LineWidth',1.5);
%             plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    yoffset=yoffset+maxy;
end
set(h1,'Visible','on');
% set(h2,'Visible','on');

% figure(); 
% imagesc(forimscB);
% colormap(gray);
% title('Blue Current');
% figure();
% imagesc(forimscR);
% colormap(gray);
% title('Red Next');

end


function out=getPvaluePSTHforEachUnit(data,prefStimset,sigStim,x1,arrayVec)

binsize=7;
plotSingleUnitData=1;
out=[];
colorSigCutoff=0.05;
% colorSigCutoff=0.0001;
plotSingleTrialData=1;

psth1vs2=nan(length(data),size(downSampMatrix(data{1},binsize),2));
psth1vs3=nan(length(data),size(downSampMatrix(data{1},binsize),2));
psth2vs3=nan(length(data),size(downSampMatrix(data{1},binsize),2));
psth1vsall=nan(length(data),size(downSampMatrix(data{1},binsize),2));
psth2vsall=nan(length(data),size(downSampMatrix(data{1},binsize),2));
psth3vsall=nan(length(data),size(downSampMatrix(data{1},binsize),2));
subx=downSampMatrix(x1(x1>=0 & x1<=1.5),binsize);
for i=1:length(data)
    d=data{i};
    sub1=downSampMatrix(d(:,x1>=0 & x1<=1.5),binsize);
    sub2=downSampMatrix(d(:,x1>=1.5 & x1<=3),binsize);
    sub3=downSampMatrix(d(:,x1>=3 & x1<=4.5),binsize);
    for j=1:length(subx)
        a=sub1(:,j);
        b=sub2(:,j);
        c=sub3(:,j);
        if any(a(~isnan(a))~=0)==0 && any(b(~isnan(b))~=0)==0 
            p=1;
        else
            p=ranksum(a(~isnan(a)),b(~isnan(b)));
        end
        if nanmean(a)>nanmean(b)
            psth1vs2(i,j)=p;
        else
            psth1vs2(i,j)=-p;
        end
        if any(a(~isnan(a))~=0)==0 && any(b(~isnan(b))~=0)==0 && any(c(~isnan(c))~=0)==0
            p=1;
        else
            p=ranksum(a(~isnan(a)),[b(~isnan(b)); c(~isnan(c))]); % vs all
        end
        if nanmean(a(~isnan(a)))>nanmean([b(~isnan(b)); c(~isnan(c))])
            psth1vsall(i,j)=p;
        else
            psth1vsall(i,j)=-p;
        end
        if any(a(~isnan(a))~=0)==0 && any(c(~isnan(c))~=0)==0
            p=1;
        else
            p=ranksum(a(~isnan(a)),c(~isnan(c)));
        end
        if nanmean(a(~isnan(a)))>nanmean(c(~isnan(c)))
            psth1vs3(i,j)=p;
        else
            psth1vs3(i,j)=-p;
        end
        if any(a(~isnan(a))~=0)==0 && any(b(~isnan(b))~=0)==0 && any(c(~isnan(c))~=0)==0
            p=1;
        else
            p=ranksum(b(~isnan(b)),[a(~isnan(a)); c(~isnan(c))]); % vs all
        end
        if nanmean(b(~isnan(b)))>nanmean([a(~isnan(a)); c(~isnan(c))])
            psth2vsall(i,j)=p;
        else
            psth2vsall(i,j)=-p;
        end
        if any(b(~isnan(b))~=0)==0 && any(c(~isnan(c))~=0)==0
            p=1;
        else
            p=ranksum(b(~isnan(b)),c(~isnan(c)));
        end
        if nanmean(b(~isnan(b)))>nanmean(c(~isnan(c)))
            psth2vs3(i,j)=p;
        else
            psth2vs3(i,j)=-p;
        end
        if any(a(~isnan(a))~=0)==0 && any(b(~isnan(b))~=0)==0 && any(c(~isnan(c))~=0)==0
            p=1;
        else
            p=ranksum(c(~isnan(c)),[a(~isnan(a)); b(~isnan(b))]); % vs all
        end
        if nanmean(c(~isnan(c)))>nanmean([a(~isnan(a)); b(~isnan(b))])
            psth3vsall(i,j)=p;
        else
            psth3vsall(i,j)=-p;
        end
    end
end

cmap1=[0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1];
% cmap2=[1      0.9   0.8  0.7  0.6 0.5 0.4 0.3 0.2 0.1];
cmap2=linspace(1,0.2,10);
% cmap2=[1      0.95   0.9  0.85  0.8 0.75 0.7 0.65 0.6 0.55];
% Make current vs next figure
h1=figure(1); 
hax=axes();
set(gcf,'Visible','off');
title('Current (blue) vs Next (red)');
h2=figure(2);
hax2=axes();
set(gcf,'Visible','off');
yoffset=0;
xoffset=subx(1)-(subx(2)-subx(1));
subx=[xoffset subx];
forimscR=zeros(length(data),length(subx)-1);
forimscB=zeros(length(data),length(subx)-1);
for i=1:length(data)
    xoffset=subx(1)-(subx(2)-subx(1));
    d=data{i};
    sub1=nanmean(downSampMatrix(d(:,x1>=0 & x1<=1.5),binsize),1);
    sub2=nanmean(downSampMatrix(d(:,x1>=1.5 & x1<=3),binsize),1);
    sub3=nanmean(downSampMatrix(d(:,x1>=3 & x1<=4.5),binsize),1);
    maxy=max([sub1 sub2 sub3],[],2);
    % Stim 1 vs 2
    if arrayVec(i,1)>=0 & 1-arrayVec(i,1)<=colorSigCutoff % Fires more for stim 1 than 2, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,1)<=0 & 1+arrayVec(i,1)<=colorSigCutoff % Fires more for stim 2 than 1, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth1vsall(i,j)=psth1vs2(i,j);
        if psth1vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 1
            cmask=[0 0 1];
        elseif psth1vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 1 
            cmask=[0 0 1];
        elseif psth1vsall(i,j)>=0 & currSig==1 & upward==0 % If FR increase and increase means stim 2
            cmask=[1 0 0];
        elseif psth1vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 2
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        colorsatp=abs(psth1vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if j==1
            figure(1); 
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub1(j) sub1(j)],'Color',currc,'LineWidth',1.5);
            hold on;
            figure(2);
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            hold on;
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else 
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub1(j-1) sub1(j)],'Color',currc,'LineWidth',1.5); 
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    % Stim 2 vs 3
    xoffset=1.5;
    if arrayVec(i,3)>=0 & 1-arrayVec(i,3)<=colorSigCutoff % Fires more for stim 2 than 3, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,3)<=0 & 1+arrayVec(i,3)<=colorSigCutoff % Fires more for stim 3 than 2, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth2vsall(i,j)=psth2vs3(i,j);
        if psth2vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 2
            cmask=[0 0 1];
        elseif psth2vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 2 
            cmask=[0 0 1];
        elseif psth2vsall(i,j)>=0 & currSig==1 & upward==0 % If FR increase and increase means stim 3
            cmask=[1 0 0];
        elseif psth2vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 3
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        colorsatp=abs(psth2vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if j==1
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub2(j) sub2(j)],'Color',currc,'LineWidth',1.5);
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub2(j-1) sub2(j)],'Color',currc,'LineWidth',1.5); 
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    % Stim 1 vs 3 (or 3 vs 1, correct order)
    xoffset=3;
    if arrayVec(i,2)>=0 & 1-arrayVec(i,2)<=colorSigCutoff % Fires more for stim 1 than 3, is significant
        upward=1;
        currSig=1;
    elseif arrayVec(i,2)<=0 & 1+arrayVec(i,2)<=colorSigCutoff % Fires more for stim 3 than 1, is significant
        upward=0;
        currSig=1;
    else
        upward=1;
        currSig=0;
    end
    for j=1:length(subx)-1
%         psth3vsall(i,j)=-psth1vs3(i,j);
        if psth3vsall(i,j)>=0 & currSig==1 & upward==0 % If FR relative increase during 3 and increase means stim 3
            cmask=[0 0 1];
        elseif psth3vsall(i,j)<=0 & currSig==1 & upward==1 % If FR decrease and decrease means stim 3 
            cmask=[0 0 1];
        elseif psth3vsall(i,j)>=0 & currSig==1 & upward==1 % If FR increase and increase means stim 1
            cmask=[1 0 0];
        elseif psth3vsall(i,j)<=0 & currSig==1 & upward==0 % If FR decrease and decrease means stim 1
            cmask=[1 0 0];
        else % Not significant change consistent with unit coding scheme
            cmask=[0 0 0];
        end
        colorsatp=abs(psth3vsall(i,j));
        cind=find(cmap1<=colorsatp,1,'last');
        colorsat=cmap2(cind);
        currc=(ones(1,3).*colorsat).*cmask;
        if j==1
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub3(j) sub3(j)],'Color',currc,'LineWidth',1.5);
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        else
            plot(hax,xoffset+[subx(j) subx(j+1)],yoffset+[sub3(j-1) sub3(j)],'Color',currc,'LineWidth',1.5);
            plot(hax2,xoffset+[subx(j) subx(j+1)],[yoffset yoffset],'Color',currc,'LineWidth',5);
            forimscR(i,j)=forimscR(i,j)+currc(1);
            forimscB(i,j)=forimscB(i,j)+currc(3);
        end
    end
    yoffset=yoffset+maxy;
end
set(h1,'Visible','on');
set(h2,'Visible','on');

figure(); 
imagesc(forimscB);
colormap(gray);
title('Blue Current');
figure();
imagesc(forimscR);
colormap(gray);
title('Red Next');

end

function out=bootstrapSum(data,prefStimset,sigStim,x1,arrayVec)

% Transform each trial into a summed representation
% weightedTrials1=addUnitsByTuning(data,prefStimset,1,sigStim,arrayVec);
% weightedTrials2=addUnitsByTuning(data,prefStimset,2,sigStim,arrayVec);
% weightedTrials3=addUnitsByTuning(data,prefStimset,3,sigStim,arrayVec);
weightedTrials1=addUnitsByTuning(data,prefStimset,1,sigStim,[]);
weightedTrials2=addUnitsByTuning(data,prefStimset,2,sigStim,[]);
weightedTrials3=addUnitsByTuning(data,prefStimset,3,sigStim,[]);

% Bootstrap weighted representation across trials
nRuns=100; % number of runs of bootstrap
takeFraction=0.25;
weightedRep1=nan(nRuns,size(data{1},2));
weightedRep2=nan(nRuns,size(data{1},2));
weightedRep3=nan(nRuns,size(data{1},2));
takeN=floor(size(weightedTrials1,1)*takeFraction);
for i=1:nRuns
    if mod(i,1000)==0
        disp(i);
    end
    currTake=1+floor(rand(1,takeN).*(size(weightedTrials1,1)-1));
    weightedRep1(i,:)=nanmean(weightedTrials1(currTake,:),1);
    weightedRep2(i,:)=nanmean(weightedTrials2(currTake,:),1);
    weightedRep3(i,:)=nanmean(weightedTrials3(currTake,:),1);
end
% Norm weightedRep for each stim cond
norm1=nanmean(weightedRep1,1);
weightedRep1=weightedRep1-nanmin(norm1);
norm1=norm1-nanmin(norm1);
weightedRep1=weightedRep1./nanmax(norm1);
norm1=norm1./nanmax(norm1);
weightedRep1=weightedRep1-nanmean(norm1);
norm2=nanmean(weightedRep2,1);
weightedRep2=weightedRep2-nanmin(norm2);
norm2=norm2-nanmin(norm2);
weightedRep2=weightedRep2./nanmax(norm2);
norm2=norm2./nanmax(norm2);
weightedRep2=weightedRep2-nanmean(norm2);
norm3=nanmean(weightedRep3,1);
weightedRep3=weightedRep3-nanmin(norm3);
norm3=norm3-nanmin(norm3);
weightedRep3=weightedRep3./nanmax(norm3);
norm3=norm3./nanmax(norm3);
weightedRep3=weightedRep3-nanmean(norm3);


% Organize into current, next and previous
currTrials=weightedRep1(:,x1>=0 & x1<=1.5);
currTrials=[currTrials; weightedRep2(:,x1>=1.5 & x1<=3)];
currTrials=[currTrials; weightedRep3(:,x1>=3 & x1<=4.5)];
nextTrials=weightedRep3(:,x1>=0 & x1<=1.5);
nextTrials=[nextTrials; weightedRep1(:,x1>=1.5 & x1<=3)];
nextTrials=[nextTrials; weightedRep2(:,x1>=3 & x1<=4.5)];
prevTrials=weightedRep2(:,x1>=0 & x1<=1.5);
prevTrials=[prevTrials; weightedRep3(:,x1>=1.5 & x1<=3)];
prevTrials=[prevTrials; weightedRep1(:,x1>=3 & x1<=4.5)];

d=5;
newcurr=downSampMatrix(currTrials,d);
newnext=downSampMatrix(nextTrials,d);
newprev=downSampMatrix(prevTrials,d);
newx=downSampAv(x1(x1>=0 & x1<=1.5),d);
out.currTrials=currTrials;
out.nextTrials=nextTrials;
out.prevTrials=prevTrials;
out.x=x1(x1>=0 & x1<=1.5);

% Plot results of boostrap
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newnext,'r',hax);
title('Sum Bootstrap');
% hl=plotLineAndStd(newx,newprev,'g',hax);
end

function aligned=meanAlignChunks(data)



end

function trialSum=addUnitsByTuning(data,prefStimset,currs,sigStim,arrayVec)

% subData=data(prefStimset==currs & sigStim<10^-5);
% subData=data(prefStimset==currs & sigStim<0.01);
% prc=prctile(sigStim,9);

if isempty(arrayVec)
    % GOOD
    prc=prctile(sigStim,9);
    % prc=prctile(sigStim(prefStimset==currs),9);
    disp(prc);
    subData=data(prefStimset==currs & sigStim<prc);
    trialSum=zeros(size(data{1}));
    for i=1:length(subData)
        trialSum=trialSum+subData{i};
    end
else
    % TRY NEW
    % Transform each trial into a weighted representation
    x=linspace(0,4.5,size(data{1},2));
%     p_cutoff=0.05;
    p_cutoff=prctile(sigStim,10);
    disp(p_cutoff);
    weightedTrials1=zeros(size(data{1}));
    weightedTrials2=zeros(size(data{1}));
    weightedTrials3=zeros(size(data{1}));
    for i=1:length(data)
        c=arrayVec(i,1);
        p1=nan(size(c));
        p1(c>=0)=1-c(c>=0);
        p1(c<0)=1+c(c<0);
        c=arrayVec(i,2);
        p2=nan(size(c));
        p2(c>=0)=1-c(c>=0);
        p2(c<0)=1+c(c<0);
        c=arrayVec(i,3);
        p3=nan(size(c));
        p3(c>=0)=1-c(c>=0);
        p3(c<0)=1+c(c<0);
        allp=[p1 p2 p3];
        for j=1:length(allp)
            if allp(j)>=p_cutoff
                continue
            else
                cForMax=arrayVec(i,j);
                d=data{i};
                if j==1
                    if cForMax>=0 
                        weightedTrials1(:,x>=0 & x<=3)=weightedTrials1(:,x>=0 & x<=3)+d(:,x>=0 & x<=3);
                    else
                        weightedTrials2(:,x>=0 & x<=3)=weightedTrials2(:,x>=0 & x<=3)+d(:,x>=0 & x<=3);
                    end
                elseif j==2
                    if cForMax>=0
                        weightedTrials1(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5))=weightedTrials1(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5))+d(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5));
                    else
                        weightedTrials3(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5))=weightedTrials3(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5))+d(:,(x>=0 & x<=1.5) | (x>=3 & x<=4.5));
                    end
                elseif j==3
                    if cForMax>=0
                        weightedTrials2(:,x>=1.5 & x<=4.5)=weightedTrials2(:,x>=1.5 & x<=4.5)+d(:,x>=1.5 & x<=4.5);
                    else
                        weightedTrials3(:,x>=1.5 & x<=4.5)=weightedTrials3(:,x>=1.5 & x<=4.5)+d(:,x>=1.5 & x<=4.5);
                    end
                end
            end
        end
    end
    if currs==1
        trialSum=weightedTrials1;
    elseif currs==2
        trialSum=weightedTrials2;
    elseif currs==3
        trialSum=weightedTrials3;
    end
end

end

function out=bootstrapWeighted(data,arrayVec,x1)

% Transform each trial into a weighted representation
weightedTrials1=zeros(size(data{1}));
weightedTrials2=zeros(size(data{1}));
weightedTrials3=zeros(size(data{1}));
for i=1:length(data)
    c=arrayVec(i,1);
    weightedTrials1=weightedTrials1+data{i}.*c;
    weightedTrials2=weightedTrials2-data{i}.*c;
    c=arrayVec(i,2);
    weightedTrials1=weightedTrials1+data{i}.*c;
    weightedTrials3=weightedTrials3-data{i}.*c;
    c=arrayVec(i,3);
    weightedTrials2=weightedTrials2+data{i}.*c;
    weightedTrials3=weightedTrials3-data{i}.*c;
%     c=arrayVec(i,1); % MESSED UP
%     weightedTrials3=weightedTrials3+data{i}.*c;
%     weightedTrials2=weightedTrials2-data{i}.*c;
%     c=arrayVec(i,2);
%     weightedTrials3=weightedTrials3+data{i}.*c;
%     weightedTrials1=weightedTrials1-data{i}.*c;
%     c=arrayVec(i,3);
%     weightedTrials2=weightedTrials2+data{i}.*c;
%     weightedTrials1=weightedTrials1-data{i}.*c;
end

% Bootstrap weighted representation across trials
nRuns=100; % number of runs of bootstrap
takeFraction=0.25;
weightedRep1=nan(nRuns,size(data{1},2));
weightedRep2=nan(nRuns,size(data{1},2));
weightedRep3=nan(nRuns,size(data{1},2));
takeN=floor(size(weightedTrials1,1)*takeFraction);
for i=1:nRuns
    if mod(i,1000)==0
        disp(i);
    end
    currTake=1+floor(rand(1,takeN).*(size(weightedTrials1,1)-1));
    weightedRep1(i,:)=nanmean(weightedTrials1(currTake,:),1);
    weightedRep2(i,:)=nanmean(weightedTrials2(currTake,:),1);
    weightedRep3(i,:)=nanmean(weightedTrials3(currTake,:),1);
end
% Norm weightedRep for each stim cond
norm1=nanmean(weightedRep1,1);
weightedRep1=weightedRep1-nanmin(norm1);
norm1=norm1-nanmin(norm1);
weightedRep1=weightedRep1./nanmax(norm1);
norm1=norm1./nanmax(norm1);
weightedRep1=weightedRep1-nanmean(norm1);
norm2=nanmean(weightedRep2,1);
weightedRep2=weightedRep2-nanmin(norm2);
norm2=norm2-nanmin(norm2);
weightedRep2=weightedRep2./nanmax(norm2);
norm2=norm2./nanmax(norm2);
weightedRep2=weightedRep2-nanmean(norm2);
norm3=nanmean(weightedRep3,1);
weightedRep3=weightedRep3-nanmin(norm3);
norm3=norm3-nanmin(norm3);
weightedRep3=weightedRep3./nanmax(norm3);
norm3=norm3./nanmax(norm3);
weightedRep3=weightedRep3-nanmean(norm3);

% Organize into current, next and previous
currTrials=weightedRep1(:,x1>=0 & x1<=1.5);
currTrials=[currTrials; weightedRep2(:,x1>=1.5 & x1<=3)];
currTrials=[currTrials; weightedRep3(:,x1>=3 & x1<=4.5)];
nextTrials=weightedRep3(:,x1>=0 & x1<=1.5);
nextTrials=[nextTrials; weightedRep1(:,x1>=1.5 & x1<=3)];
nextTrials=[nextTrials; weightedRep2(:,x1>=3 & x1<=4.5)];
prevTrials=weightedRep2(:,x1>=0 & x1<=1.5);
prevTrials=[prevTrials; weightedRep3(:,x1>=1.5 & x1<=3)];
prevTrials=[prevTrials; weightedRep1(:,x1>=3 & x1<=4.5)];

d=5;
newcurr=downSampMatrix(currTrials,d);
newnext=downSampMatrix(nextTrials,d);
newprev=downSampMatrix(prevTrials,d);
newx=downSampAv(x1(x1>=0 & x1<=1.5),d);
out.currTrials=currTrials;
out.nextTrials=nextTrials;
out.prevTrials=prevTrials;
out.x=x1(x1>=0 & x1<=1.5);

% Plot results of boostrap
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newnext,'r',hax);
% hl=plotLineAndStd(newx,newprev,'g',hax);
end

function hl=plotLineAndErr(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1)./sqrt(size(data,1)),'y',hax,hl);

end

function hl=plotLineAndStd(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1),'y',hax,hl);

end


function temp=makeArraySameLengthAsVec(a,b)

% Makes vec b same length as vec a

temp=nan(1,length(a));
if length(a)>size(b,2)
    temp(:,1:size(b,2))=b;
elseif size(b,2)>length(a)
    temp=b(:,1:length(a));
else
    temp=b;
end

end