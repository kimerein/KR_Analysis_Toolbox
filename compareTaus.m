function [taus]=compareTaus(xpoints,datas)

showTestFigs=0;
showExpFits=0;
smoothCurves=1;
nSmooths=1;
LEDonsets=0;
topNormWindows=[-0.02 0];
bottomNormWindows=[0.08 0.1];
w_default=0.003;
tau_default=0.01;
colorOrder={'k','r','b','g','m','c',[0.1 0.1 0.1],[0.3 0.3 0.3],[0.5 0.5 0.5]};
% colorOrder={'k','r','b','g','m','c'};
smoothWindow=10;

% Normalize all shut-off curves
shutoffDatas=cell(length(datas),1);
for i=1:length(datas)
    allCurves=datas{i};
    topNormWindow=topNormWindows;
    bottomNormWindow=bottomNormWindows;
    subXpoints{i}=xpoints(xpoints>=topNormWindow(1) & xpoints<=bottomNormWindow(2));
    currRunningData=zeros(size(allCurves,1),length(subXpoints{i}));
    for j=1:size(allCurves,1) 
        subCurve=allCurves(j,xpoints>=topNormWindow(1) & xpoints<=bottomNormWindow(2));
        subCurve=subCurve-mean(subCurve(subXpoints{i}>=bottomNormWindow(1) & subXpoints{i}<=bottomNormWindow(2)));
        topVal=mean(subCurve(subXpoints{i}>=topNormWindow(1) & subXpoints{i}<=topNormWindow(2)));
        subCurve=subCurve./topVal;
        if smoothCurves==1
            for k=1:nSmooths
                subCurve=smooth(subCurve,smoothWindow);
            end
        end
        currRunningData(j,:)=subCurve;
    end
    figure();
    for j=1:size(allCurves,1)
        if j>length(colorOrder)
            plot(subXpoints{1},currRunningData(j,:),'Color',colorOrder{1});
        else
            plot(subXpoints{1},currRunningData(j,:),'Color',colorOrder{j});
        end
        hold on;
    end
%     figure(); 
%     plot(subXpoints{1},currRunningData');
    shutoffDatas{i}=currRunningData;
end

% Get non-parametric data on shut-offs
% There is noise in data, so get both the index of first crossing on way
% down and the index of first threshold crossing on the way up from zero in
% negative time -- take the average of these two times as the threshold
% crossing time
% Top is 1, bottom is 0
thresh=0.75;
percent75times=cell(length(datas),1);
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    currCrossingData=zeros(size(allCurves,1),1);
    for j=1:size(allCurves,1)
        fixedCurve=allCurves(j,:);
        fixedCurve(1:(find(subXpoints{i}>LEDonsets,1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets,1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets;
        if showTestFigs==1
            figure(); 
            plot(s,allCurves(j,:),'Color','k');
            hold on; 
            line([s(downCross) s(downCross)],[0 1]);
            line([s(upCross) s(upCross)],[0 1]);
            line([crossTime crossTime],[0 1]);
        end
%         disp(i);
%         disp(crossTime);
%         disp(LEDonsets);
    end
    percent75times{i}=currCrossingData;
end
thresh=0.5;
percent50times=cell(length(datas),1);
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    currCrossingData=zeros(size(allCurves,1),1);
    for j=1:size(allCurves,1)
        fixedCurve=allCurves(j,:);
        fixedCurve(1:(find(subXpoints{i}>LEDonsets,1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets,1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets;
    end
    percent50times{i}=currCrossingData;
end     
thresh=0.25;
percent25times=cell(length(datas),1);
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    currCrossingData=zeros(size(allCurves,1),1);
    for j=1:size(allCurves,1)
        fixedCurve=allCurves(j,:);
        fixedCurve(1:(find(subXpoints{i}>LEDonsets,1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets,1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets;
    end
    percent25times{i}=currCrossingData;
end

% Get exponential fits and parametric summaries of shut-offs
ws=cell(length(datas),1);
taus=cell(length(datas),1);
if showTestFigs==1
    showExpFits=1;
end
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    topNormWindow=topNormWindows;
    bottomNormWindow=bottomNormWindows;
    currWsData=zeros(size(allCurves,1),1);
    currTausData=zeros(size(allCurves,1),1);
    for j=1:size(allCurves,1)
        [fTau,fW,fR,fi]=fitExponential_forParametric(subXpoints{i},allCurves(j,:),LEDonsets,0,0.003,0.01,bottomNormWindow,showExpFits);
        if showExpFits==1
            choice=questdlg('Accept fit?', 'Yes', 'No');
            switch choice
                case 'Yes'
                    currWsData(j)=fW;
                    close(fi);
                case 'No'
                    currWsData(j)=-10000;
                    close(fi);
                case 'Cancel'
                    return
            end
        else
            currWsData(j)=fW;
        end
        [fTau,fW,fR,fi]=fitExponential_forParametric(subXpoints{i},allCurves(j,:),LEDonsets,1,0.003,0.01,bottomNormWindow,showExpFits);
        if showExpFits==1
            choice=questdlg('Accept fit?', 'Yes', 'No');
            switch choice
                case 'Yes'
                    currTausData(j)=fTau;
                    close(fi);
                case 'No'
                    currTausData(j)=-10000;
                    close(fi);
                case 'Cancel'
                    return
            end
        else
            currTausData(j)=fTau;
        end
    end
    ws{i}=currWsData;
    taus{i}=currTausData;
end

% Plot results
figure(); 
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    plot(subXpoints{1},mean(allCurves,1));
    hold all;
end

% Plot summary metrics
figure(); 
subplot(1,5,1);
for i=1:length(percent75times)
    alls=percent75times{i};
    for j=1:size(alls,1)
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('ms to 0.75','FontSize',8);
%     hold on;
end
subplot(1,5,2);
for i=1:length(percent50times)
    alls=percent50times{i};
    for j=1:size(alls,1)
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('ms to 0.5','FontSize',8);
%     hold on;
end
subplot(1,5,3);
for i=1:length(percent25times)
    alls=percent25times{i};
    for j=1:size(alls,1)
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('ms to 0.25','FontSize',8);
%     hold on;
end
subplot(1,5,4);
for i=1:length(ws)
    alls=ws{i};
    if showExpFits==1
        alls=alls(alls~=-10000);
        c=colorOrder(alls~=-10000);
    else
        c=colorOrder;
    end
    for j=1:size(alls,1)
%         scatter(i,alls(j)*1000,[],colorOrder{j});
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('w (ms)','FontSize',8);
%     hold on;
end
subplot(1,5,5);
for i=1:length(taus)
    alls=taus{i};
    if showExpFits==1
        alls=alls(alls~=-10000);
        c=colorOrder(alls~=-10000);
    else
        c=colorOrder;
    end
    for j=1:size(alls,1)
%         scatter(i,alls(j)*1000,[],colorOrder{j});
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('tau (ms)','FontSize',8);
%     hold on;
end
