function compareTauMetrics(xpoints,datas,LEDonsets)
% LEDonsets should be cell array
showTestFigs=0;
smoothCurves=1;
nSmooths=2;
% LEDonset=1.3;
% topNormWindows={[1.18 1.2]; [1.18 1.2]};
% topNormWindows={[1.25 1.3]; [1.25 1.3]};
% topNormWindows={[1.33 1.34]; [1.33 1.34]};
% topNormWindows={[1.27 1.34]; [1.27 1.34]};
% bottomNormWindows={[1.43 1.44]; [1.43 1.44]};
% bottomNormWindows={[1.4 1.45]; [1.4 1.45]};
% topNormWindows={[3.16 3.2]; [3.26 3.3]};
% bottomNormWindows={[3.28 3.32]; [3.38 3.42]};
topNormWindows={[1.28 1.3]; [1.28 1.3]};
bottomNormWindows={[1.36 1.38]; [1.36 1.38]};
w_default=0.003;
% w_default=0;
tau_default=0.01;
colorOrder={'k','r','b','g','m','c','f',[0.1 0.1 0.1],[0.3 0.3 0.3],[0.5 0.5 0.5]};
smoothWindow=12;

% Normalize all shut-off curves
shutoffDatas=cell(length(datas),1);
for i=1:length(datas)
    allCurves=datas{i};
    topNormWindow=topNormWindows{i};
    bottomNormWindow=bottomNormWindows{i};
    subXpoints{i}=xpoints(xpoints>=topNormWindow(1) & xpoints<=bottomNormWindow(2));
    currRunningData=zeros(size(allCurves,1),length(subXpoints{i}));
%     topNormWindow=topNormWindows{1};
%     bottomNormWindow=bottomNormWindows{1};
%     subXpoints{i}=xpoints(xpoints>=topNormWindow(1) & xpoints<=bottomNormWindow(2));
%     currRunningData=zeros(size(allCurves,1),length(subXpoints{i}));
    for j=1:size(allCurves,1) 
%         topNormWindow=topNormWindows{j};
%         bottomNormWindow=bottomNormWindows{j};
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
    plot(subXpoints{1},currRunningData');
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
        fixedCurve(1:(find(subXpoints{i}>LEDonsets{i},1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets{i},1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets{i};
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
%         disp(LEDonsets{i});
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
        fixedCurve(1:(find(subXpoints{i}>LEDonsets{i},1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets{i},1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets{i};
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
        fixedCurve(1:(find(subXpoints{i}>LEDonsets{i},1)-1))=ones(length(1:(find(subXpoints{i}>LEDonsets{i},1)-1)),1);
        downCross=find(fixedCurve<thresh,1);
        upCross=find(fixedCurve(end:-1:1)>thresh,1);
        upCross=length(fixedCurve)-(upCross-1);
%         downCross=find(allCurves(j,:)<thresh,1);
%         upCross=find(allCurves(j,end:-1:1)>thresh,1);
%         upCross=length(allCurves(j,:))-(upCross-1);
        avCross=ceil(mean([downCross upCross]));
        s=subXpoints{i};
        crossTime=s(avCross);
        currCrossingData(j)=crossTime-LEDonsets{i};
    end
    percent25times{i}=currCrossingData;
end

% Get exponential fits and parametric summaries of shut-offs
ws=cell(length(datas),1);
taus=cell(length(datas),1);
for i=1:length(shutoffDatas)
    allCurves=shutoffDatas{i};
    topNormWindow=topNormWindows{i};
    bottomNormWindow=bottomNormWindows{i};
    currWsData=zeros(size(allCurves,1),1);
    currTausData=zeros(size(allCurves,1),1);
    for j=1:size(allCurves,1)
        [fTau,fW,fR]=fitExponential_forParametric(subXpoints{i},allCurves(j,:),LEDonsets{i},0,0.003,0.01,bottomNormWindow,showTestFigs);
        currWsData(j)=fW;
        [fTau,fW,fR]=fitExponential_forParametric(subXpoints{i},allCurves(j,:),LEDonsets{i},1,0.003,0.01,bottomNormWindow,showTestFigs);
        currTausData(j)=fTau;
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
    disp(alls);
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
    for j=1:size(alls,1)
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
    disp(alls);
    for j=1:size(alls,1)
        scatter(i,alls(j)*1000,[],colorOrder{j});
        hold on;
    end
%     scatter(ones(size(alls,1),1)*i,alls*1000,[],colorOrder{i});
    xlabel('tau conds','FontSize',8);
    ylabel('tau (ms)','FontSize',8);
%     hold on;
end
