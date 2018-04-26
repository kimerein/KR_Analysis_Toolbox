function putTogetherAlphaPhaseOutput(acrossExpts)

bootout.currTrials=[];
bootout.nextTrials=[];
bootout.prevTrials=[];
bootout.x=[];
outsum.currTrials=[];
outsum.nextTrials=[];
outsum.prevTrials=[];
outsum.x=[];
currTrials=[];
nextTrials=[];
prevTrials=[];
bins=[];
p=[];
p_prev=[];
currSum=[];
nextSum=[];
prevSum=[];
for i=1:length(acrossExpts)
    a=acrossExpts{i};
    bootout.currTrials=[bootout.currTrials; a.bootout.currTrials];
    bootout.nextTrials=[bootout.nextTrials; a.bootout.nextTrials];
    bootout.prevTrials=[bootout.prevTrials; a.bootout.prevTrials];
    bootout.x=[bootout.x; a.bootout.x];
    outsum.currTrials=[outsum.currTrials; a.outsum.currTrials];
    outsum.nextTrials=[outsum.nextTrials; a.outsum.nextTrials];
    outsum.prevTrials=[outsum.prevTrials; a.outsum.prevTrials];
    outsum.x=[outsum.x; a.outsum.x];
    currTrials=[currTrials; a.currTrials];
    nextTrials=[nextTrials; a.nextTrials];
    prevTrials=[prevTrials; a.prevTrials];
    currSum(i,:)=nansum(a.currTrials,1);
    nextSum(i,:)=nansum(a.nextTrials,1);
    prevSum(i,:)=nansum(a.prevTrials,1);
    bins=[bins; a.outbins];
    p=[p; a.outps];
    p_prev=[p_prev; a.outps_prev];
end

% Plot average summed significantly increased units
x=bootout.x(1,:);
currSum=currSum-repmat(min(currSum(:,x>=1 & x<=1.5),[],2),1,size(currSum,2));
currSum=currSum./repmat(max(currSum(:,x>=1 & x<=1.5),[],2),1,size(currSum,2));
nextSum=nextSum-repmat(min(nextSum(:,x>=1 & x<=1.5),[],2),1,size(nextSum,2));
nextSum=nextSum./repmat(max(nextSum(:,x>=1 & x<=1.5),[],2),1,size(nextSum,2));
prevSum=prevSum-repmat(min(prevSum(:,x>=1 & x<=1.5),[],2),1,size(prevSum,2));
prevSum=prevSum./repmat(max(prevSum(:,x>=1 & x<=1.5),[],2),1,size(prevSum,2));
d=2;
newcurr1=downSampMatrix(currSum,d);
newnext1=downSampMatrix(nextSum,d);
newprev1=downSampMatrix(prevSum,d);
newx1=downSampAv(bootout.x(1,:),d);
figure(); 
hax=axes();
% plot(newx1,newcurr1','Color','k');
% hold on; 
% plot(newx1,newnext1','Color','r');
hl=plotLineAndErr(newx1,newcurr1,'k',hax);
hold on; 
hl=plotLineAndErr(newx1,newnext1,'r',hax);
title('Raw Sum Output: Curr vs Next');
figure(); 
hax=axes();
% plot(newx1,newcurr1','Color','k');
% hold on; 
% plot(newx1,newprev1','Color','g');
hl=plotLineAndErr(newx1,newcurr1,'k',hax);
hold on; 
hl=plotLineAndErr(newx1,newprev1,'g',hax);
title('Raw Sum Output: Curr vs Prev');

% Plot output of bootstrap
d=5;
newcurr=downSampMatrix(bootout.currTrials,d);
newnext=downSampMatrix(bootout.nextTrials,d);
newprev=downSampMatrix(bootout.prevTrials,d);
newx=downSampAv(bootout.x(1,:),d);
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newnext,'r',hax);
title('Bootstrap Output: Curr vs Next');
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newprev,'g',hax);
title('Bootstrap Output: Curr vs Prev');

% Plot output of sum bootstrap
d=5;
newcurr=downSampMatrix(outsum.currTrials,d);
newnext=downSampMatrix(outsum.nextTrials,d);
newprev=downSampMatrix(outsum.prevTrials,d);
newx=downSampAv(outsum.x(1,:),d);
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newnext,'r',hax);
title('Bootstrap Sum: Curr vs Next');
figure(); 
hax=axes();
hl=plotLineAndErr(newx,newcurr,'k',hax);
hold on; 
hl=plotLineAndErr(newx,newprev,'g',hax);
title('Bootstrap Sum: Curr vs Prev');

% Plot example trials
d=2;
newCurrTrials=downSampMatrix(currTrials,d);
newNextTrials=downSampMatrix(nextTrials,d);
newPrevTrials=downSampMatrix(prevTrials,d);
nx=downSampAv(bootout.x(1,:),d);
figure(); 
plot(nx(nx>=1.0 & nx<=1.5),newNextTrials(:,nx>=1.0 & nx<=1.5)'+rand(size(newNextTrials(:,nx>=1.0 & nx<=1.5)))'.*20,'Color','r');
hold on;
plot(nx(nx>=1.0 & nx<=1.5),newCurrTrials(:,nx>=1.0 & nx<=1.5)'+rand(size(newCurrTrials(:,nx>=1.0 & nx<=1.5)))'.*20,'Color','k');
title('Example Trials Curr vs Next');
figure(); 
plot(nx(nx>=1.0 & nx<=1.5),newPrevTrials(:,nx>=1.0 & nx<=1.5)'+rand(size(newPrevTrials(:,nx>=1.0 & nx<=1.5)))'.*20,'Color','g');
hold on;
plot(nx(nx>=1.0 & nx<=1.5),newCurrTrials(:,nx>=1.0 & nx<=1.5)'+rand(size(newCurrTrials(:,nx>=1.0 & nx<=1.5)))'.*20,'Color','k');
title('Example Trials Curr vs Prev');

% Plot sum across all trials
figure(); 
plot(nx(nx>=1.0 & nx<=1.5),nansum(newNextTrials(:,nx>=1.0 & nx<=1.5),1),'Color','r');
hold on;
plot(nx(nx>=1.0 & nx<=1.5),nansum(newCurrTrials(:,nx>=1.0 & nx<=1.5),1),'Color','k');
title('Sum Curr vs Next');
figure(); 
plot(nx(nx>=1.0 & nx<=1.5),nansum(newPrevTrials(:,nx>=1.0 & nx<=1.5),1),'Color','g');
hold on;
plot(nx(nx>=1.0 & nx<=1.5),nansum(newCurrTrials(:,nx>=1.0 & nx<=1.5),1),'Color','k');
title('Sum Curr vs Prev');

% Plot pvals
figure();
for i=1:size(p,1)
    scatter(bins(i,:),p(i,:));
    hold all;
end
title('Pvals Curr vs Next');
figure();
for i=1:size(p,1)
    scatter(bins(i,:),p_prev(i,:));
    hold all;
end
title('Pvals Curr vs Prev');

% Recalculate pvals using data from all experiments
d=1;
newcurr=downSampMatrix(currTrials,d);
newnext=downSampMatrix(nextTrials,d);
newprev=downSampMatrix(prevTrials,d);
newx=downSampAv(x,d);
offsets=[0 0.05 0.1];
% offsets=[0 0.025 0.05 0.075];
tryBins=sort(1.5:-0.15:0);
% tryBins=sort(1.5:-0.1:0);
allp=nan(length(offsets),length(tryBins)-1);
o=cell(1,length(offsets));
allp_prev=nan(length(offsets),length(tryBins)-1);
for j=1:length(offsets)
    tryBins=sort(1.5-offsets(j):-0.2:0);
    o{j}=tryBins;
    pvalsAcrossTrial=ones(1,length(tryBins)-1);
    pvalsAcrossTrial_prev=ones(1,length(tryBins)-1);
    for i=1:length(tryBins)-1
        window=[tryBins(i) tryBins(i+1)];
        pvalsAcrossTrial(i)=compareDistributions(newx,newnext,newcurr,window);
        pvalsAcrossTrial_prev(i)=compareDistributions(newx,newprev,newcurr,window);
    end
    allp(j,end-length(pvalsAcrossTrial)+1:end)=pvalsAcrossTrial;
    allp_prev(j,end-length(pvalsAcrossTrial_prev)+1:end)=pvalsAcrossTrial_prev;
end
figure(); 
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
title('Pvals Across Expts');
disp(allp);
figure(); 
scatter(outbins,outps_prev);
title('Pvals Prev Across Expts');

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

function hl=plotLineAndErr(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1)./sqrt(size(data,1)),'y',hax,hl);

end

function hl=plotLineAndStd(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1),'y',hax,hl);

end
    