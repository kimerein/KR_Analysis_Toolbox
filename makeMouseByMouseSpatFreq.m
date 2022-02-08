function makeMouseByMouseSpatFreq(mouseName,t)

%load('Z:\Kim\FF_manuscript\revisions J neuro\Fig 3 spatial freqs\raw F1 PSTHs\mouseA_point03sf_noThetaNoLED.mat')

doNorm=false;

sfnames={'_point02sf', '_point03sf', '_point04sf', '_point08sf'};
% sf=[0.02, 0.03, 0.04, 0.08];
sf=[0.02, 0.03, 0.04, 0.05]; % just for plotting equally spaced
% sfnames={'_point02sf', '_point04sf', '_point08sf'};
% sf=[0.02, 0.04, 0.08];
% sfnames={'_point02sf', '_point03sf', '_point04sf', '_point06sf'};

allsfev=nan(length(sfnames),4);
unitbyunit=nan(length(sfnames),4,1000);
unitbyunitspont=nan(length(sfnames),4,1000);
backupMouseName=mouseName;
baseWindow=[3 3.75];
% baseWindow=[3.5 3.75];
for i=1:length(sfnames)
    sfname=sfnames{i};
    if ~isempty(regexp(sfname,'_point03sf')) & ~isempty(regexp(mouseName,'mouseF'))
        mouseName=[mouseName(1:end-1) 'A'];  
    elseif ~isempty(regexp(sfname,'_point03sf')) & ~isempty(regexp(mouseName,'mouseC'))
        mouseName=[mouseName(1:end-1) 'A'];
    else
        mouseName=backupMouseName;
    end
    if ~isempty(regexp(sfname,'_point08sf')) & ~isempty(regexp(mouseName,'mouseE'))
        sfname='_point06sf';
    end
    a=load([mouseName sfname '_noThetaNoLED.mat']);
    temp=a.temp;
    evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t>baseWindow(1) & t<=baseWindow(2)),2);
    spont=nanmean(temp(:,t>3.25 & t<=3.75),2);
    noThetaNoLED_spont=spont;
    noThetaNoLED_evoked=evoked;
    a=load([mouseName sfname '_thetaNoLED.mat']);
    temp=a.temp;
    evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t>baseWindow(1) & t<=baseWindow(2)),2);
    spont=nanmean(temp(:,t>3.25 & t<=3.75),2);
    thetaNoLED_spont=spont;
    thetaNoLED_evoked=evoked;
    if ~isempty(regexp(sfname,'_point03sf')) & ~isempty(regexp(mouseName,'mouseE'))
%         sc=12; 
        sc=3;
        noThetaNoLED_evoked=noThetaNoLED_evoked*sc;
        thetaNoLED_evoked=thetaNoLED_evoked*sc;
    end
    p=signrank(thetaNoLED_evoked-noThetaNoLED_evoked);
    [f,x,flo,fup]=ecdf(thetaNoLED_evoked-noThetaNoLED_evoked,'Bounds','on');
    figure();
    plot(x,f,'Color','k');
    hold on;
    disp([sfname 'theta no LED vs no theta no LED, pval ' num2str(p)]);
    title(sfname);
    
    a=load([mouseName sfname '_noThetaLED.mat']);
    temp=a.temp;
    evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t>baseWindow(1) & t<=baseWindow(2)),2);
    spont=nanmean(temp(:,t>3.25 & t<=3.75),2);
    noThetaLED_spont=spont;
    noThetaLED_evoked=evoked;
    if ~isempty(regexp(sfname,'_point03sf')) & ~isempty(regexp(mouseName,'mouseE'))
        noThetaLED_evoked=noThetaLED_evoked*sc;
    end
    p=signrank(noThetaLED_evoked-noThetaNoLED_evoked);
    [f,x,flo,fup]=ecdf(noThetaLED_evoked-noThetaNoLED_evoked,'Bounds','on');
    plot(x,f,'Color','c');
    disp([sfname 'no theta no LED vs no theta LED, pval ' num2str(p)]);
    
    a=load([mouseName sfname '_thetaLED.mat']);
    temp=a.temp;
    evoked=nanmean(temp(:,t>4 & t<=6.5),2)-nanmean(temp(:,t>baseWindow(1) & t<=baseWindow(2)),2);
    spont=nanmean(temp(:,t>3.25 & t<=3.75),2);
    thetaLED_spont=spont;
    thetaLED_evoked=evoked;
    if ~isempty(regexp(sfname,'_point03sf')) & ~isempty(regexp(mouseName,'mouseE'))
        thetaLED_evoked=thetaLED_evoked*sc;
    end
    
    if ~isempty(regexp(sfname,'_point08sf')) & ~isempty(regexp(mouseName,'mouseD'))
        thetaLED_evoked([65 70 71])=nan;
    end
    p=signrank(thetaLED_evoked-thetaNoLED_evoked);
    [f,x,flo,fup]=ecdf(thetaLED_evoked-thetaNoLED_evoked,'Bounds','on');
    plot(x,f,'Color','b');
    disp([sfname 'theta no LED vs theta LED, pval ' num2str(p)]);
    
    p=signrank(thetaLED_evoked-noThetaLED_evoked);
    [f,x,flo,fup]=ecdf(thetaLED_evoked-noThetaLED_evoked,'Bounds','on');
    plot(x,f,'Color','m');
    disp([sfname 'theta LED vs no theta LED, pval ' num2str(p)]);

%     visuallyresponsive=noThetaNoLED_evoked>0 | thetaNoLED_evoked>0 | noThetaLED_evoked>0 | thetaLED_evoked>0; 
%     allsfev(i,1)=nanmean(noThetaNoLED_evoked(visuallyresponsive==1));
%     allsfev(i,2)=nanmean(thetaNoLED_evoked(visuallyresponsive==1));
%     allsfev(i,3)=nanmean(noThetaLED_evoked(visuallyresponsive==1));
%     allsfev(i,4)=nanmean(thetaLED_evoked(visuallyresponsive==1));
    allsfev(i,1)=nanmean(noThetaNoLED_evoked(noThetaNoLED_evoked>0));
    allsfev(i,2)=nanmean(thetaNoLED_evoked(thetaNoLED_evoked>0));
    allsfev(i,3)=nanmean(noThetaLED_evoked(noThetaLED_evoked>0));
    allsfev(i,4)=nanmean(thetaLED_evoked(thetaLED_evoked>0));
    unitbyunit(i,1,1:length(noThetaNoLED_evoked(noThetaNoLED_evoked>0)))=noThetaNoLED_evoked(noThetaNoLED_evoked>0);
    unitbyunit(i,2,1:length(thetaNoLED_evoked(thetaNoLED_evoked>0)))=thetaNoLED_evoked(thetaNoLED_evoked>0);
    unitbyunit(i,3,1:length(noThetaLED_evoked(noThetaLED_evoked>0)))=noThetaLED_evoked(noThetaLED_evoked>0);
    unitbyunit(i,4,1:length(thetaLED_evoked(thetaLED_evoked>0)))=thetaLED_evoked(thetaLED_evoked>0);
    unitbyunitspont(i,1,1:length(noThetaNoLED_spont))=noThetaNoLED_spont;
    unitbyunitspont(i,2,1:length(thetaNoLED_spont))=thetaNoLED_spont;
    unitbyunitspont(i,3,1:length(noThetaLED_spont))=noThetaLED_spont;
    unitbyunitspont(i,4,1:length(thetaLED_spont))=thetaLED_spont;
end

if doNorm
    normVals=nanmean(unitbyunit(:,1,:),1);
    normed_unitbyunit=unitbyunit./repmat(normVals(1:end),size(unitbyunit,1),size(unitbyunit,2),1);
    figure();
    plot(sf,nanmean(normed_unitbyunit(:,1,:),3),'Color',[0.5 0.5 0.5]);
    hold on;
    plot(sf,nanmean(normed_unitbyunit(:,2,:),3),'Color','k');
    plot(sf,nanmean(normed_unitbyunit(:,3,:),3),'Color','c');
    plot(sf,nanmean(normed_unitbyunit(:,4,:),3),'Color','b');
else
    figure();
    plot(sf,allsfev(:,1),'Color',[0.5 0.5 0.5]);
    hold on;
    plotStderrs(sf,allsfev(:,1),unitbyunit,1,[0.5 0.5 0.5],'false');
    plot(sf,allsfev(:,2),'Color','k');
    plotStderrs(sf,allsfev(:,2),unitbyunit,2,'k','false');
    plot(sf,allsfev(:,3),'Color','c');
    plotStderrs(sf,allsfev(:,3),unitbyunit,3,'c','false');
    plot(sf,allsfev(:,4),'Color','b');
    plotStderrs(sf,allsfev(:,4),unitbyunit,4,'b','false');
    unitbyunitspont(1,:,:)=nanmean(unitbyunitspont,1);
    scatter(0.1,nanmean(unitbyunitspont(1,1,:),3),[],[0.5 0.5 0.5]);
    plotStderrs(0.1,nanmean(unitbyunitspont(1,1,:),3),unitbyunitspont,1,[0.5 0.5 0.5],true);
    scatter(0.1,nanmean(unitbyunitspont(1,2,:),3),[],'k');
    plotStderrs(0.1,nanmean(unitbyunitspont(1,2,:),3),unitbyunitspont,1,'k',true);
    scatter(0.1,nanmean(unitbyunitspont(1,3,:),3),[],'c');
    plotStderrs(0.1,nanmean(unitbyunitspont(1,3,:),3),unitbyunitspont,1,'c',true);
    scatter(0.1,nanmean(unitbyunitspont(1,4,:),3),[],'b');
    plotStderrs(0.1,nanmean(unitbyunitspont(1,4,:),3),unitbyunitspont,1,'b',true);
end
xlim([0.01 0.11]);

end

function plotStderrs(sf,me,unitbyunit,whichind,c,justonce)

for i=1:size(unitbyunit,1)
    line([sf(i) sf(i)],[me(i)-nanstd(unitbyunit(i,whichind,:),[],3)./nansum(~isnan(unitbyunit(i,whichind,:))) me(i)+nanstd(unitbyunit(i,whichind,:),[],3)./nansum(~isnan(unitbyunit(i,whichind,:)))],'Color',c);
    if justonce==true
        break
    end
end

end