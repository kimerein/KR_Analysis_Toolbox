function [fig_handles,axes_handles,colorMapCode]=plotLFPconds(data,Fs,stimsForSweeps,ledForSweeps,params,averageOrAdd,colorMapCode,matchNumTrials,showSignificance,statsBinWidth,showPhoto,useLED,photoSignal,ledSignal)


fig_handles=[];
axes_handles=[];
colorMapCode=[];

if any(size(data)==0)
    disp('No data!');
    return
end

if floor(statsBinWidth/2)<=0
    binOffset=1;
else
    binOffset=floor(statsBinWidth/2);
end

f1=figure;
fig_handles=[fig_handles; f1];

% If no LED data, assume LED is off for all trials
if isempty(ledForSweeps)
    ledForSweeps=zeros(size(data,1),1);
end
if isempty(stimsForSweeps)
    stimsForSweeps=ones(size(data,1),1);
end

% Plan figure layout
nfigRows=3;
nfigCols=1;

% If colorMapCode does not exist, make one
% Assign different colors to different stimulus conditions
if isempty(colorMapCode)
    [lineColors,colorMapCode]=makeColorCode(params);
end

% Normalize ledForSweeps to make it a logical array
% Assumes a step LED (that is, only 2 LED values)
minLED=min(ledForSweeps);
maxLED=max(ledForSweeps);
if minLED==maxLED
    ledForSweeps(ledForSweeps==0)=0;
    ledForSweeps(ledForSweeps~=0)=1;
else
    ledForSweeps(ledForSweeps==minLED)=0;
    ledForSweeps(ledForSweeps==maxLED)=1;
end

% Plot the signal collapsed over all stimulus conditions
curr_h=axesmatrix(nfigRows,nfigCols,1);
% set1=logical(ledForSweeps);
% set2=~logical(ledForSweeps);
set1b=find(logical(ledForSweeps)==1);
set2b=find(~logical(ledForSweeps)==1);
if matchNumTrials==1 && useLED
    [set1b,set2b]=matchNumberOfTrials(set1b,set2b);
end
set1=ismember(1:length(ledForSweeps),set1b);
set2=ismember(1:length(ledForSweeps),set2b);
disp('length(set1)');
disp(length(set1b));
disp('length(set2)');
disp(length(set2b));
times=0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength;
if any(set1~=0)
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color','b');
end
hold on;
if any(set2~=0)
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color','k');
end
if ~any(set1~=0)
    minall=min(combineData(averageOrAdd,data(set2,400:end-400),1));
    maxall=max(combineData(averageOrAdd,data(set2,400:end-400),1));
elseif ~any(set2~=0)
    minall=min(combineData(averageOrAdd,data(set1,400:end-400),1));
    maxall=max(combineData(averageOrAdd,data(set1,400:end-400),1));
else
    min1=min(combineData(averageOrAdd,data(set1,400:end-400),1));
    max1=max(combineData(averageOrAdd,data(set1,400:end-400),1));
    min2=min(combineData(averageOrAdd,data(set2,400:end-400),1));
    max2=max(combineData(averageOrAdd,data(set2,400:end-400),1));
    minall=min(min1,min2);
    maxall=max(max1,max2);
end
sigSignal_lessThanPoint05=[];
sigSignal_lessThanPoint01=[];
sigSignal_lessThanPoint001=[];
if showSignificance==1
    [sigSignal_lessThanPoint05,sigSignal_lessThanPoint01,sigSignal_lessThanPoint001,indsForBins]=ttestSignals(data(set1,:),data(set2,:),statsBinWidth,binOffset);
    minall=minall-abs(maxall-minall)/30;
end
for i=1:length(sigSignal_lessThanPoint05)
    if sigSignal_lessThanPoint001(i)~=0
        if i==length(sigSignal_lessThanPoint001)
            line([times(indsForBins(i)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.1 0.1 0.1],'LineWidth',6);
        else
            line([times(indsForBins(i)) times(indsForBins(i+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.1 0.1 0.1],'LineWidth',6);
        end
        continue
    end
    if sigSignal_lessThanPoint01(i)~=0
        if i==length(sigSignal_lessThanPoint01)
            line([times(indsForBins(i)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.5 0.5 0.5],'LineWidth',6);
        else
            line([times(indsForBins(i)) times(indsForBins(i+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.5 0.5 0.5],'LineWidth',6);
        end
        continue
    end
    if sigSignal_lessThanPoint05(i)~=0
        if i==length(sigSignal_lessThanPoint05)
            line([times(indsForBins(i)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.75 0.75 0.75],'LineWidth',6);
        else
            line([times(indsForBins(i)) times(indsForBins(i+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.75 0.75 0.75],'LineWidth',6);
        end
    end
end
title('Av. Stim.-Triggered Signal');
s=struct('figure_handle',curr_h,'figure_description','total signal collapsed over all stimulus conditions');
axes_handles=[axes_handles; s];
sigRange=abs(maxall-minall);
if showPhoto==1
    photoRange=max(photoSignal)-min(photoSignal);
    a=(1/20)*(sigRange/photoRange);
    photoSignal2=a*photoSignal-min(a*photoSignal);
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,photoSignal2+(minall-(1/9.5)*sigRange),'Color',[0.2 0.2 0.2]);
    minall=minall-(1/9.5)*sigRange;
end
%useLED==0;
if useLED==1
    ledRange=max(ledSignal)-min(ledSignal);
    a=(1/20)*(sigRange/ledRange);
    ledSignal2=a*ledSignal-min(a*ledSignal);
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,ledSignal2+(minall-(1/9.5)*sigRange),'Color','blue');
    minall=minall-(1/9.5)*sigRange;
end
axis([0 params.totalTrialLength minall maxall]);

% Plot the signal collapsed over Stimulus Variable 1
curr_h=axesmatrix(nfigRows,nfigCols,2);
trialsForConds=cell(length(params.Var2_values),1);
for i=1:length(params.Var2_values)
    theseConds=i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values);
    trialsForConds{i}=find(ismember(stimsForSweeps,theseConds));
end
trialsForConds=matchManyNumberOfTrials(trialsForConds);
disp('trialsForConds');
disp(trialsForConds);
for i=1:length(params.Var2_values)
    %theseConds=i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values);
    cOffset=length(params.Var1_values)*length(params.Var2_values)+length(params.Var1_values)+i;
%     set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
%     set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
%     if matchNumTrials==1
%         [set1,set2]=matchNumberOfTrials(set1,set2);
%     end
%     plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color',colorMapCode(cOffset).color);
%     hold on;
%     plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color','k');
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(trialsForConds{i},:),1),'Color',colorMapCode(cOffset).color);
    hold on;
end
title(['Stim.-Triggered Signal - Collapsed over ' params.Var1_name]);
s=struct('figure_handle',curr_h,'figure_description',['Superimposed: total signal collapsed over ' params.Var1_name]);
axis tight;
axes_handles=[axes_handles; s];

% Plot the signal collapsed over Stimulus Variable 2
curr_h=axesmatrix(nfigRows,nfigCols,3);
trialsForConds=cell(length(params.Var1_values),1);
for i=1:length(params.Var1_values)
    theseConds=i:1:i+length(params.Var2_values)-1;
    trialsForConds{i}=find(ismember(stimsForSweeps,theseConds));
end
trialsForConds=matchManyNumberOfTrials(trialsForConds);
disp('trialsForConds');
disp(trialsForConds);
for i=1:length(params.Var1_values)
    %theseConds=i:1:i+length(params.Var2_values)-1;
    cOffset=length(params.Var1_values)*length(params.Var2_values)+i;
%     set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
%     set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
%     if matchNumTrials==1
%         [set1,set2]=matchNumberOfTrials(set1,set2);
%     end
%     plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color',colorMapCode(cOffset).color);
%     hold on;
%     plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color','k');
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(trialsForConds{i},:),1),'Color',colorMapCode(cOffset).color);
    hold on;
end
title(['Stim.-Triggered Signal - Collapsed over ' params.Var2_name]);
s=struct('figure_handle',curr_h,'figure_description',['Superimposed: total signal collapsed over ' params.Var2_name]);
axis tight;
axes_handles=[axes_handles; s];


% Now separate plots for all combinations of stim. Var1 and stim. Var2
% conditions
f2=figure;
fig_handles=[fig_handles; f2];
% Plan figure layout
nfigRows=1+length(params.Var1_values);
nfigCols=1+length(params.Var2_values);

k=2;
% Collapsed over Var1
keepTheseAxes=[];
ccs=cell(2,1);
for i=1:length(params.Var2_values) 
    curr_h=axesmatrix(nfigRows,nfigCols,k);
    keepTheseAxes=[keepTheseAxes; curr_h];
    theseConds=i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values);
    cOffset=length(params.Var1_values)*length(params.Var2_values)+length(params.Var1_values)+i;
    set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
    set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
    set1b=find(set1==1);
    set2b=find(set2==1);
    if matchNumTrials==1 && useLED==1
        [set1b,set2b]=matchNumberOfTrials(set1b,set2b);
    end
    set1=ismember(1:length(ledForSweeps),set1b);
    set2=ismember(1:length(ledForSweeps),set2b);
    ccs{1}='k';
    ccs{2}='k';
    if useLED==1
        ccs{1}=colorMapCode(cOffset).color;
        ccs{2}='k';
    else
        ccs{1}='k';
        ccs{2}=colorMapCode(cOffset).color;
    end
    if ~isempty(set1)
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color',ccs{1});
    end
    hold on;
    if ~isempty(set2)
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color',ccs{2});
    end
    k=k+1;
    s=struct('figure_handle',curr_h,'figure_description',['Separate: total signal collapsed over ' params.Var1_name]);
    axes_handles=[axes_handles; s];
    if i==1
        if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            maxSignal=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            minSignal=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
        else
            maxSignal=0;
            minSignal=0;
        end
        if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
            maxSignal2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            minSignal2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
        else
            maxSignal2=maxSignal;
            minSignal2=minSignal;
        end
        maxSignal=max(maxSignal,maxSignal2);
        minSignal=min(minSignal,minSignal2);
        max2=maxSignal; 
        max_2=maxSignal;
        min2=minSignal;
        min_2=minSignal;
    else
        if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            max2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            min2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
        end
        if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
            max_2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            min_2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
        end
        if ~any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            max2=max_2;
            min2=min_2;
        elseif ~any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
        else
            max2=max(max2,max_2);
            min2=min(min2,min_2);
        end
    end
    sigSignal_lessThanPoint05=[];
    sigSignal_lessThanPoint01=[];
    sigSignal_lessThanPoint001=[];
    sigRange=abs(max2-min2);
    if showSignificance==1
        [sigSignal_lessThanPoint05,sigSignal_lessThanPoint01,sigSignal_lessThanPoint001,indsForBins]=ttestSignals(data(set1,:),data(set2,:),statsBinWidth,binOffset);
        min2=min2-(abs(max2-min2))/50;
    end
    for j=1:length(sigSignal_lessThanPoint05)
        if sigSignal_lessThanPoint001(j)~=0
            if j==length(sigSignal_lessThanPoint001)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
            end
            continue
        end
        if sigSignal_lessThanPoint01(j)~=0
            if j==length(sigSignal_lessThanPoint01)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
            end
            continue
        end
        if sigSignal_lessThanPoint05(j)~=0
            if j==length(sigSignal_lessThanPoint05)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.7 0.7 0.7],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.7 0.7 0.7],'LineWidth',4);
            end
        end
    end
    if showPhoto==1
        photoRange=max(photoSignal)-min(photoSignal);
        a=(1/20)*(sigRange/photoRange);
        photoSignal2=a*photoSignal-min(a*photoSignal);
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,photoSignal2+(min2-(1/9.5)*sigRange),'Color',[0.2 0.2 0.2]);
        min2=min2-(1/9.5)*sigRange;
    end
    if useLED==1
        ledRange=max(ledSignal)-min(ledSignal);
        a=(1/20)*(sigRange/ledRange);
        ledSignal2=a*ledSignal-min(a*ledSignal);
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,ledSignal2+(min2-(1/9.5)*sigRange),'Color','blue');
        min2=min2-(1/9.5)*sigRange;
    end
    if max2>maxSignal
        maxSignal=max2;
    end
    if min2<minSignal
        minSignal=min2;
    end
end
if useLED
    for i=1:length(keepTheseAxes) 
        axis(keepTheseAxes(i),[0 params.totalTrialLength minSignal maxSignal]);
    end
end

% Collapsed over Var2
keepTheseAxes=[];
ccs=cell(2,1);
for i=1:length(params.Var1_values)
    curr_h=axesmatrix(nfigRows,nfigCols,i*nfigCols+1);
    keepTheseAxes=[keepTheseAxes; curr_h];
    theseConds=i:1:i+length(params.Var2_values)-1;
    cOffset=length(params.Var1_values)*length(params.Var2_values)+i;
    set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
    set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
    set1b=find(set1==1);
    set2b=find(set2==1);
    if matchNumTrials==1 && useLED==1
        [set1b,set2b]=matchNumberOfTrials(set1b,set2b);
    end
    set1=ismember(1:length(ledForSweeps),set1b);
    set2=ismember(1:length(ledForSweeps),set2b);
    ccs{1}='k';
    ccs{2}='k';
    if useLED==1
        ccs{1}=colorMapCode(cOffset).color;
        ccs{2}='k';
    else
        ccs{1}='k';
        ccs{2}=colorMapCode(cOffset).color;
    end
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color',ccs{1});
    hold on;
    plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color',ccs{2});
    s=struct('figure_handle',curr_h,'figure_description',['Separate: total signal collapsed over ' params.Var2_name]);
    axes_handles=[axes_handles; s];
    if i==1
        if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            maxSignal=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            minSignal=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
        else
            maxSignal=0;
            minSignal=0;
        end
        if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
            maxSignal2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            minSignal2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
        else
            maxSignal2=maxSignal;
            minSignal2=minSignal;
        end
        maxSignal=max(maxSignal,maxSignal2);
        minSignal=min(minSignal,minSignal2);
        max2=maxSignal; 
        max_2=maxSignal;
        min2=minSignal;
        min_2=minSignal;
    else
        if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            max2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            min2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
        end
        if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
            max_2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            min_2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
        end
        if ~any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
            max2=max_2;
            min2=min_2;
        elseif ~any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
        else
            max2=max(max2,max_2);
            min2=min(min2,min_2);
        end
    end
    sigSignal_lessThanPoint05=[];
    sigSignal_lessThanPoint01=[];
    sigSignal_lessThanPoint001=[];
    sigRange=abs(max2-min2);
    if showSignificance==1
        [sigSignal_lessThanPoint05,sigSignal_lessThanPoint01,sigSignal_lessThanPoint001,indsForBins]=ttestSignals(data(set1,:),data(set2,:),statsBinWidth,binOffset);
        min2=min2-(abs(max2-min2))/50;
    end
    for j=1:length(sigSignal_lessThanPoint05)
        if sigSignal_lessThanPoint001(j)~=0
            if j==length(sigSignal_lessThanPoint001)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
            end
            continue
        end
        if sigSignal_lessThanPoint01(j)~=0
            if j==length(sigSignal_lessThanPoint01)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
            end
            continue
        end
        if sigSignal_lessThanPoint05(j)~=0
            if j==length(sigSignal_lessThanPoint05)
                line([times(indsForBins(j)) times(end)],[min2 min2],'Color',[0.75 0.75 0.75],'LineWidth',4);
            else
                line([times(indsForBins(j)) times(indsForBins(j+1))],[min2 min2],'Color',[0.75 0.75 0.75],'LineWidth',4);
            end
        end
    end
    if showPhoto==1
        photoRange=max(photoSignal)-min(photoSignal);
        a=(1/20)*(sigRange/photoRange);
        photoSignal2=a*photoSignal-min(a*photoSignal);
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,photoSignal2+(min2-(1/9.5)*sigRange),'Color',[0.2 0.2 0.2]);
        min2=min2-(1/9.5)*sigRange;
    end
    if useLED==1
        ledRange=max(ledSignal)-min(ledSignal);
        a=(1/20)*(sigRange/ledRange);
        ledSignal2=a*ledSignal-min(a*ledSignal);
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,ledSignal2+(min2-(1/9.5)*sigRange),'Color','blue');
        min2=min2-(1/9.5)*sigRange;
    end
    if max2>maxSignal
        maxSignal=max2;
    end
    if min2<minSignal
        minSignal=min2;
    end
end
if useLED
    for i=1:length(keepTheseAxes)
        axis(keepTheseAxes(i),[0 params.totalTrialLength minSignal maxSignal]);
    end
end

% Separate plots for each stimulus condition (combo. of Var1 and Var2)
k=nfigCols+2;
keepTheseAxes=[];
ccs=cell(2,1);
for i=1:length(params.Var1_values)
    for j=1:length(params.Var2_values)
        curr_h=axesmatrix(nfigRows,nfigCols,k);
        keepTheseAxes=[keepTheseAxes curr_h];
        theseConds=(i-1)*length(params.Var2_values)+j;
        cOffset=(i-1)*length(params.Var2_values)+j;
        set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
        set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
        set1b=find(set1==1);
        set2b=find(set2==1);
        if matchNumTrials==1 && useLED==1
            [set1b,set2b]=matchNumberOfTrials(set1b,set2b);
        end
        set1=ismember(1:length(ledForSweeps),set1b);
        set2=ismember(1:length(ledForSweeps),set2b);
        ccs{1}='k';
        ccs{2}='k';
        if useLED==1
            ccs{1}=colorMapCode(cOffset).color;
            ccs{2}='k';
        else
            ccs{1}='k';
            ccs{2}=colorMapCode(cOffset).color;
        end
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set1,:),1),'Color',ccs{1});
        hold on;
        plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,combineData(averageOrAdd,data(set2,:),1),'Color',ccs{2});
        title([num2str(params.Var1_values(i)) ' ' num2str(params.Var2_values(j))]);
        s=struct('figure_handle',curr_h,'figure_description',['Average: ' num2str(params.Var1_values(i)) ' ' num2str(params.Var2_values(j))]);
        axes_handles=[axes_handles; s];
        if j==1 && i==1
            if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
                maxSignal=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
                minSignal=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            else
                maxSignal=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds),400:end-400),1));
                minSignal=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds),400:end-400),1));
            end
            if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
                maxSignal2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
                minSignal2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            else
                maxSignal2=maxSignal;
                minSignal2=minSignal;
            end
            maxSignal=max(maxSignal,maxSignal2);
            minSignal=min(minSignal,minSignal2);
            max2=maxSignal;
            max_2=maxSignal;
            min2=minSignal;
            min_2=minSignal;
        else
            if any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
                max2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
                min2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&ledForSweeps,400:end-400),1));
            end
            if any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
                max_2=max(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
                min_2=min(combineData(averageOrAdd,data(ismember(stimsForSweeps,theseConds)&~ledForSweeps,400:end-400),1));
            end
            if ~any(ismember(stimsForSweeps,theseConds)&ledForSweeps)
                max2=max_2;
                min2=min_2;
            elseif ~any(ismember(stimsForSweeps,theseConds)&~ledForSweeps)
            else
                max2=max(max2,max_2);
                min2=min(min2,min_2);
            end
        end
        sigSignal_lessThanPoint05=[];
        sigSignal_lessThanPoint01=[];
        sigSignal_lessThanPoint001=[];
        sigRange=abs(max2-min2);
        if showSignificance==1
            [sigSignal_lessThanPoint05,sigSignal_lessThanPoint01,sigSignal_lessThanPoint001,indsForBins]=ttestSignals(data(set1,:),data(set2,:),statsBinWidth,binOffset);
            min2=min2-(abs(max2-min2))/50;
        end
        for q=1:length(sigSignal_lessThanPoint05)
            if sigSignal_lessThanPoint001(q)~=0
                if q==length(sigSignal_lessThanPoint001)
                    line([times(indsForBins(q)) times(end)],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
                else
                    line([times(indsForBins(q)) times(indsForBins(q+1))],[min2 min2],'Color',[0.1 0.1 0.1],'LineWidth',4);
                end
                continue
            end
            if sigSignal_lessThanPoint01(q)~=0
                if q==length(sigSignal_lessThanPoint01)
                    line([times(indsForBins(q)) times(end)],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
                else
                    line([times(indsForBins(q)) times(indsForBins(q+1))],[min2 min2],'Color',[0.5 0.5 0.5],'LineWidth',4);
                end
                continue
            end
            if sigSignal_lessThanPoint05(q)~=0
                if q==length(sigSignal_lessThanPoint05)
                    line([times(indsForBins(q)) times(end)],[min2 min2],'Color',[0.75 0.75 0.75],'LineWidth',4);
                else
                    line([times(indsForBins(q)) times(indsForBins(q+1))],[min2 min2],'Color',[0.75 0.75 0.75],'LineWidth',4);
                end
            end
        end
        if showPhoto==1
            photoRange=max(photoSignal)-min(photoSignal); 
            a=(1/20)*(sigRange/photoRange);
            photoSignal2=a*photoSignal-min(a*photoSignal);
            plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,photoSignal2+(min2-(1/9.5)*sigRange),'Color',[0.2 0.2 0.2]);
            min2=min2-(1/9.5)*sigRange;
        end
        if useLED==1    
            ledRange=max(ledSignal)-min(ledSignal);
            a=(1/20)*(sigRange/ledRange);
            ledSignal2=a*ledSignal-min(a*ledSignal);
            plot(0:params.totalTrialLength/(size(data,2)-1):params.totalTrialLength,ledSignal2+(min2-(1/9.5)*sigRange),'Color','blue');
            min2=min2-(1/9.5)*sigRange;
        end
        if mod(k,nfigCols)==0
            k=k+2;
        else
            k=k+1;
        end 
        if max2>maxSignal
            maxSignal=max2;
        end
        if min2<minSignal
            minSignal=min2;
        end
    end
end
if useLED
    for i=1:length(keepTheseAxes)
        axis(keepTheseAxes(i),[0 params.totalTrialLength minSignal maxSignal]);
    end
end
end

function r=combineData(averageOrAdd,data,dim)
    if strcmp(averageOrAdd,'average')
        r=mean(data,dim);
    else
        r=sum(data,dim);
    end
end