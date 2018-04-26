function [togetherLayers,normLayers,unitShutoffs,allCSDs,allHalfMaxTimes]=putLayerDataTogether(useDir,listing)
%listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');

subtractBaseForUnits=1;
weightAverage=1;
plotThalDepAndInd=0;
nHistBins=3;
colorByDepth=0;
excludeEverywhere=1;
exclude=[]; %15
% useWeights=[1 2 3 3 3 3 3];
useWeights=[];
% exclude=[8 16 8-2 6-2 9-2 11-2 14-2 15-2 17-2]; % Set A
% exclude=[8 16 8-2 3-2 4-2 5-2 7-2 10-2 12-2 13-2 16-2 18-2]; % Set B
doSpont=0; % let zero be zero
timeBefore=0.5; % time before LED onset to take
timeAfter=0.5; % time after LED onset to take
% baseRelativeToLEDonset=[-1 -0.5]; 
% baseRelativeToLEDonset=[-0.2 -0.05]; 
% baseRelativeToLEDonset=[0 1]; 
baseRelativeToLEDonset=[-0.5 -0.25];
% peakRelativeToLEDonset=[-0.3 0];
peakRelativeToLEDonset=[-0.15 0];
% downSampFactor=30;
downSampFactor=1;
downSampRaw=1;
alpha1=0.05;
alpha2=0.01;

cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};

% listing=dir(useDir);
% Put PSTH data together
indBefore=timeBefore*1000;
indAfter=timeAfter*1000;
useBase=baseRelativeToLEDonset.*1000;
usePeak=peakRelativeToLEDonset.*1000;
% useBase=baseRelativeToLEDonset.*1000+timeBefore*1000;
% usePeak=peakRelativeToLEDonset.*1000+timeBefore*1000;
ncon=zeros(1,length(listing));
nled=zeros(1,length(listing));
allcon=zeros(length(listing),indBefore+indAfter);
allled=zeros(length(listing),indBefore+indAfter);
useX=-indBefore+1:1:indAfter; % in ms
togetherLayers.y1=zeros(45,indBefore+indAfter);
togetherLayers.y2=zeros(45,indBefore+indAfter);
gotData=0;
unitTaus=[];
unitEvChs=[];
allCSDs=zeros(500,30000);
allHalfMaxTimes=nan(45,length(listing));
for i=1:length(listing)
    currFolder=listing(i).name;
    disp(i);
%     if ~exist([useDir '\' currFolder '\exptData.mat'],'file')
%         continue
%     end
%     a=load([useDir '\' currFolder '\exptData.mat']);
    if ~exist([useDir '\' currFolder '\layerDataFIXED.mat'],'file')...
      || ~exist([useDir '\' currFolder '\CSDalignment.mat'],'file')
        continue
    end
    gotData=gotData+1;
    a=load([useDir '\' currFolder '\layerDataFIXED.mat']);
    layerData=a.layerData;
    a=load([useDir '\' currFolder '\CSDalignment.mat']);
    CSDalignment=a.CSDalignment;
    a=load([useDir '\' currFolder '\exptData.mat']);
    exptData=a.exptData;
    a=load([useDir '\' currFolder '\paramsForSilencing.mat']);
    params=a.params;
    if ~exist([useDir '\' currFolder '\units_shutoff.mat'],'file')
        allUnitsTimeCourse.taus=[];
        allUnitsTimeCourse.calibrated_evCh=[];
    else
        a=load([useDir '\' currFolder '\units_shutoff.mat']);
        allUnitsTimeCourse=a.allUnitsTimeCourse;
    end
    if ~exist([useDir '\' currFolder '\dataInt.mat'],'file')
        dataInt=[];
    else
        a=load([useDir '\' currFolder '\dataInt.mat']);
        dataInt=a.dataInt;
    end
    ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
    halfMaxTimes=nan(size(layerData.y2,1),1);
    for j=1:size(layerData.y2,1)
        l.x=layerData.x;
        l.y2=layerData.y2(j,:);
        halfMaxTimes(j)=fitSingleLayerShutoff(l,ledOnset);
    end
    firstInd=find(exptData.xpoints>=ledOnset,1);
    allcon(i,:)=exptData.ypoints1(firstInd-indBefore+1:firstInd+indAfter);
    allled(i,:)=exptData.ypoints2(firstInd-indBefore+1:firstInd+indAfter);
    ncon(i)=exptData.numtrials1;
    nled(i)=exptData.numtrials2;
    subLayerData.y1=layerData.y1(:,firstInd-indBefore+1:firstInd+indAfter);
    subLayerData.y2=layerData.y2(:,firstInd-indBefore+1:firstInd+indAfter);
%     % Find time to half-max
%     % Have to down-sample or smooth, apparently
%     smoothed_y2=subLayerData.y2;
%     for j=1:size(smoothed_y2,1)
%         smoothed_y2(j,:)=smooth(smoothed_y2(j,:),5);
%     end
%     topVal=mean(smoothed_y2(:,1:499),2);
%     bottomVal=mean(smoothed_y2(:,580:600),2);
%     threshVals=bottomVal+(topVal-bottomVal)./2;
%     halfMaxInds=nan(size(threshVals));
%     for j=1:length(halfMaxInds)
% %         if sum(subLayerData.y2(j,:)==0)>0.7*length(subLayerData.y2(j,:))
% %             continue
% %         end
%         fromTopInd=find(smoothed_y2(j,500:600)<=threshVals(j),1,'first');
%         if isempty(fromTopInd)
%             fromTopInd=100;
%         end
%         fromBottomInd=find(smoothed_y2(j,600:-1:500)>=threshVals(j),1,'first');
%         if isempty(fromBottomInd)
%             fromBottomInd=100;
%         end
%         fromBottomInd=101-fromBottomInd;
% %         halfMaxInds(j)=mean([fromTopInd fromBottomInd]);
%         halfMaxInds(j)=mean(fromTopInd);
%     end
%     halfMaxTimes=halfMaxInds*(exptData.xpoints(2)-exptData.xpoints(1));
    allHalfMaxTimes(16-CSDalignment:1:16-CSDalignment+length(halfMaxTimes)-1,i)=halfMaxTimes;
    % CSDalignment is a scalar
    useRows=16-CSDalignment:1:16-CSDalignment+size(subLayerData.y1,1)-1;
    togetherLayers.y1(useRows,:)=togetherLayers.y1(useRows,:)+subLayerData.y1;
    togetherLayers.y2(useRows,:)=togetherLayers.y2(useRows,:)+subLayerData.y2;
    allUnitsTimeCourse.taus=allUnitsTimeCourse.taus(allUnitsTimeCourse.taus<=0.5);
    allUnitsTimeCourse.calibrated_evCh=allUnitsTimeCourse.calibrated_evCh(allUnitsTimeCourse.taus<=0.5);
    unitTaus=[unitTaus; allUnitsTimeCourse.taus];
    unitEvChs=[unitEvChs; 16-CSDalignment+allUnitsTimeCourse.calibrated_evCh'];
    if ~isempty(dataInt)
        startStim=exptData.stimulusOn(1);
%         CSDx=linspace(0,params.trialDuration,size(dataInt,2));
        CSDx=0:(1/64000):(1/64000)*(size(dataInt,2)-1);
        stimOnsetInd=find(CSDx>=startStim,1,'first');
        takeInds=stimOnsetInd-14999:stimOnsetInd+15000;
%         takeInds=CSDx>=startStim-0.25 & CSDx<=startStim+0.25;
%         allCSDs(250-floor(10*CSDalignment):250-floor(10*CSDalignment)+size(dataInt,1)-1,1:sum(takeInds))=allCSDs(250-floor(10*CSDalignment):250-floor(10*CSDalignment)+size(dataInt,1)-1,1:sum(takeInds))+dataInt(:,takeInds);
        allCSDs(250-floor(10*CSDalignment):250-floor(10*CSDalignment)+size(dataInt,1)-1,:)=allCSDs(250-floor(10*CSDalignment):250-floor(10*CSDalignment)+size(dataInt,1)-1,:)+dataInt(:,takeInds);
    end
end
togetherLayers.y1=togetherLayers.y1/gotData;
togetherLayers.y2=togetherLayers.y2/gotData;
unitShutoffs.unitTaus=unitTaus;
unitShutoffs.unitEvChs=unitEvChs;

figure(); 
takeNoNanMean=zeros(size(allHalfMaxTimes,1),1);
takeNoNanStd=zeros(size(allHalfMaxTimes,1),1);
for i=1:size(allHalfMaxTimes,1)
    takeNoNanMean(i)=mean(allHalfMaxTimes(i,~isnan(allHalfMaxTimes(i,:))));
    takeNoNanStd(i)=std(allHalfMaxTimes(i,~isnan(allHalfMaxTimes(i,:))),[],2)./sqrt(length(allHalfMaxTimes(i,~isnan(allHalfMaxTimes(i,:)))));
end
errorbar(1:length(takeNoNanMean),takeNoNanMean,takeNoNanStd);

figure(); 
imagesc(allCSDs);

figure(); 
imagesc(togetherLayers.y1);
figure(); 
imagesc(togetherLayers.y2);

figure(); 
scatter(unitEvChs,unitTaus);

normLayers.y1=togetherLayers.y1;
normLayers.y2=togetherLayers.y2;
for i=1:size(normLayers.y1,1)
    mi=min(normLayers.y1(i,:));
    normLayers.y1(i,:)=normLayers.y1(i,:)-mi;
    ma=max(normLayers.y1(i,:));
    if ma>0
        normLayers.y1(i,:)=normLayers.y1(i,:)/ma;
    end
end
for i=1:size(normLayers.y2,1)
    mi=min(normLayers.y2(i,:));
    normLayers.y2(i,:)=normLayers.y2(i,:)-mi;
    ma=max(normLayers.y2(i,:));
    if ma>0
        normLayers.y2(i,:)=normLayers.y2(i,:)/ma;
    end
end
figure(); 
imagesc(normLayers.y1);
figure(); 
imagesc(normLayers.y2);

% makeCSDtypeFig(normLayers.y1);
% makeCSDtypeFig(normLayers.y2);

% Make weighted average
% scon=zeros(1,indBefore+indAfter);
% sled=zeros(1,indBefore+indAfter);
% nusing=1;
% for i=1:length(ncon)
%     if ismember(i,exclude)
%         disp('excluding');
%         disp(listing(i).name);
%         continue
%     end
%     if weightAverage==1 && isempty(useWeights)
%         scon=scon+ncon(i)*allcon(i,:);
%         sled=sled+nled(i)*allled(i,:);
%     elseif weightAverage==1
%         scon=scon+useWeights(i)*allcon(i,:);
%         sled=sled+useWeights(i)*allled(i,:);
%     else
%         scon=scon+allcon(i,:);
%         sled=sled+allled(i,:);
%         nusing=nusing+1;
%     end
% end
% if weightAverage==1
%     scon=scon/sum(ncon);
%     sled=sled/sum(nled);
% else
%     scon=scon/nusing;
%     sled=sled/nusing;
% end
% if downSampRaw==1
%     raw.data1=downSampAv(scon,downSampFactor);
%     raw.data2=downSampAv(sled,downSampFactor);
%     raw.x=downSampAv(useX,downSampFactor)-timeBefore;
% else
%     raw.data1=scon;
%     raw.data2=sled;
%     raw.x=useX-timeBefore;
% end
% figure(); 
% plot(raw.x,raw.data1,'Color','k');
% hold on;
% plot(raw.x,raw.data2,'Color','r');
% 
% % Make normalized average
% if doSpont==1
%     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% else
%     [newX,newy1,s1,newy2,s2]=putTogetherSilencingPSTH(useX,allcon(~ismember(1:size(allcon,1),exclude),:),allled(~ismember(1:size(allcon,1),exclude),:),downSampFactor,alpha1,alpha2,doSpont,useBase,usePeak,1);
% end
% norm.x=newX;
% norm.y1=newy1;
% norm.y1_std=s1;
% norm.y2=newy2;
% norm.y2_std=s2;
% 
% % Put unit data together
% led1UnitData=[];
% led2UnitData=[];
% baseUnitData1=[];
% baseUnitData2=[];
% calibratedEvCh=[];
% waveformWidth=[];
% setOfWaveforms=zeros(1000,48);
% k=1;
% for i=1:length(listing)
%     currFolder=listing(i).name;
%     if exist([useDir '\' currFolder '\unitsFx.mat'],'file')
%         a=load([useDir '\' currFolder '\unitsFx.mat']);
%         b=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
%         if subtractBaseForUnits==1
%             c=load([useDir '\' currFolder '\baseUnitsFx.mat']);
%         end
%     else
%         continue
%     end
%     if excludeEverywhere==1 && ismember(i,exclude)
%         continue
%     end
%     unitsFx=a.unitsFx;
%     if subtractBaseForUnits==1
%         baseUnits=c.unitsFx;
%     end
%     newassignsinfo=b.newassignsinfo;
%     led1UnitData=[led1UnitData unitsFx.led1Spikes];
%     led2UnitData=[led2UnitData unitsFx.led2Spikes];
%     if subtractBaseForUnits==1
%         baseUnitData1=[baseUnitData1 baseUnits.led1Spikes];
%         baseUnitData2=[baseUnitData2 baseUnits.led2Spikes];
%     end
%     calibratedEvCh=[calibratedEvCh newassignsinfo.calibrated_evCh];
%     waveformWidth=[waveformWidth newassignsinfo.waveformWidths];
%     setOfWaveforms(k:k+size(newassignsinfo.waveforms,1)-1,:)=newassignsinfo.waveforms(:,1:48);
%     k=k+size(newassignsinfo.waveforms,1);
% end
% unitWaveforms=setOfWaveforms(1:k-1,:);
% if subtractBaseForUnits==1
%     led1UnitData=led1UnitData-baseUnitData1;
%     led2UnitData=led2UnitData-baseUnitData2;
% end
% units.led1=led1UnitData;
% units.led2=led2UnitData;
% figure();
% if colorByDepth==1
%     k=1;
%     for i=1:4:16
%         scatter(led1UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),led2UnitData(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
%         hold on;
%         k=k+1;
%     end
% else
%     scatter(led1UnitData,led2UnitData);
% end
% xlabel('LED off -- Hz');
% ylabel('LED on -- Hz');
% 
% thalDep=led1UnitData-led2UnitData;
% thalIndep=led2UnitData;
% if plotThalDepAndInd==1
%     figure();
%     if colorByDepth==1
%         k=1;
%         for i=1:4:16
%             scatter(thalDep(calibratedEvCh>=i & calibratedEvCh<=i+3),thalIndep(calibratedEvCh>=i & calibratedEvCh<=i+3),[],cs{k});
%             hold on;
%             k=k+1;
%         end
%     else
%         scatter(thalDep,thalIndep);
%     end
%     xlabel('Thal. Dependent -- Hz');
%     ylabel('Thal. Independent -- Hz');
% end
% 
% % Make histogram of fractional suppresion
% fracSupps=[];
% for i=1:length(listing)
%     if ismember(i,exclude)
%         continue
%     end
%     currFolder=listing(i).name;
%     if exist([useDir '\' currFolder '\fracSupp.mat'],'file')
%         a=load([useDir '\' currFolder '\fracSupp.mat']);
%     else
%         continue
%     end
%     fracSupp=a.fracSupp;
%     fracSupps=[fracSupps fracSupp];
% end
% [n,xout]=hist(fracSupps,nHistBins);
% figure(); 
% plot(xout,n);
% histData.xout=xout;
% histData.n=n;

end

function unitTau=fitSingleLayerShutoff(layerData,ledOnset)
% function unitTimeCourse=fitSingleLayerShutoff(spikes,useAss,exptData,defaultXLength,params)
waitS=0.003; % in s
finishS=waitS+0.06; % in s
Awindow=[-0.05 0]; % in s
zeroWindow=[0.05 0.1]; % in s

% try
%     [x,y1,y2]=scriptForComparingLayers(filtspikes(spikes,0,'assigns',useAss),exptData.AIs,[],[],exptData,params);
%     unitTimeCourse.x=x;
%     unitTimeCourse.y1=y1;
%     unitTimeCourse.y2=y2;
% catch
%     unitTimeCourse.x=NaN(1,defaultXLength);
%     unitTimeCourse.y1=NaN(1,defaultXLength);
%     unitTimeCourse.y2=NaN(1,defaultXLength);
%     unitTimeCourse.tau=nan;
%     return
% end
x=layerData.x';
y2=layerData.y2';
subX=x(x>=ledOnset+waitS & x<=ledOnset+finishS);
subX=subX-min(subX);
subY=y2(x>=ledOnset+waitS & x<=ledOnset+finishS);
forcedA=mean(y2(x>=ledOnset+Awindow(1) & x<=ledOnset+Awindow(2)));
subY=subY-mean(subY(subX>=zeroWindow(1) & subX<=zeroWindow(2)));
% Down sample if not many spikes
newSubX=subX;
newSubY=subY;
for i=2:20
    if sum(newSubY(newSubX<0.01)<=0)>0
        newSubX=downSampAv(subX,i);
        newSubY=downSampAv(subY,i);
    else
        break
    end
end
subX=newSubX;
subY=newSubY;
subY(subY==0)=0.01;
% Try fitting method 1
goodFit=0;
unitTau=nan;
try
    fo=fitoptions('Method','NonlinearLeastSquares','StartPoint',[max(subY) 0.01 max(subY)/10 0.1]);
    ft=fittype('a.*exp(-x./b)+c.*exp(-x./d)','options',fo);
    [fitout,gof,output]=fit(subX,subY,ft);
    if abs(fitout.b)>1 && fitout.d<0.1
        unitTau=fitout.d;
        goodFit=1;
    elseif abs(fitout.d)>1 && fitout.b<0.1
        unitTau=fitout.b;
        goodFit=1;
    end
catch
end
if goodFit==0
    % Try fitting method 2
    firstZero=find(subY<=0.011,1,'first');
    log_subY=log(subY(1:firstZero)); % linear if y is exponential
    log_subX=subX(1:firstZero);
    areRealEl=zeros(size(log_subY));
    for i=1:length(log_subY)
        if isreal(log_subY(i))
            areRealEl(i)=1;
        end
    end
    log_subX=log_subX(areRealEl & ~isinf(log_subY));
    log_subY=log_subY(areRealEl & ~isinf(log_subY)); 
    logForcedA=log(forcedA);
    if ~isempty(log_subY)
        fzt=mean(log_subY(2:end)-logForcedA)/(mean(log_subX(2:end))-log_subX(1));
        forcezero_tau=-1/fzt;
        forceresid=subY-subY(1)*exp(-subX/forcezero_tau);
        forceNormOfResid=(sum(forceresid.^2))^(1/2);
        % Fitting method 3
        try
            p=polyfit(log_subX(~isinf(log_subY))',log_subY(~isinf(log_subY))',1);
            tau=-1/p(1);
            A=exp(p(2));
            resid=subY-A*exp(-subX/tau);
            normOfResid=(sum(resid.^2))^(1/2);
        catch
            unitTimeCourse.tau=forcezero_tau;
            return
        end
        if normOfResid<forceNormOfResid
            unitTau=tau;
        else
            unitTau=forcezero_tau;
        end
    end
end
if unitTau<0
    unitTau=nan;
end
unitTimeCourse.tau=unitTau;
if ~isreal(unitTau)
    disp('im');
end
end

function makeCSDtypeFig(dataMatrix)
csd=zeros(size(dataMatrix,1),size(dataMatrix,2));
csd(1,:)=(2*dataMatrix(1,:)+dataMatrix(2,:))/3;
csd(end,:)=(2*dataMatrix(end,:)+dataMatrix(end-1,:))/3;
for i=2:size(dataMatrix,1)-1
    csd(i,:)=(dataMatrix(i-1,:)+2*dataMatrix(i,:)+dataMatrix(i+1,:))/4;
end
for i=1:size(dataMatrix,1)
    csd(i,:)=smooth(csd(i,:),17);
end
data1=csd(1:2:end,:);
data2=csd(2:2:end,:);
data1=diff(data1,2,1);
data2=diff(data2,2,1);
data=zeros(length(data1),2*size(data1,1));
DataCounter=1;
for i=1:size(data1,1)
    data(:,DataCounter)=data1(i,:);
    DataCounter=DataCounter+1;
    if DataCounter>size(data2,1)
        break
    end
    data(:,DataCounter)=data2(i,:);
    DataCounter=DataCounter+1;
end
data=-data';

% 2D linear interpolation
% [x y]=size(data);
% y=1:y;
% x=1:x;
% [xi yi]=meshgrid(1:1:max(x),1:0.05:max(y));
% dataInt=interp2(x,y,data',xi,yi);
[x y]=size(data);
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
dataInt=interp2(x,y,data',xi,yi);

figure(); 
dataInt=dataInt';
imagesc(dataInt);
drawnow;
end

function [newXpoints,newDataY1,stds1,newDataY2,stds2]=putTogetherSilencingPSTH(xpoints,dataY1,dataY2,downSampFactor,alpha1,alpha2,keepZero,base,peak,scale)

% Use a larger bin size for PSTHs
for i=1:size(dataY1,1)
    [newDataY1(i,:),newDataY2(i,:),newXpoints]=useLargerBinsize(dataY1(i,:),dataY2(i,:),xpoints,downSampFactor);
end

% Normalize PSTHs
baselineInds=find(newXpoints>base(1) & newXpoints<base(2));
peakInds=find(newXpoints>peak(1) & newXpoints<peak(2));
sharedBaseline=mean(mean(newDataY1(:,baselineInds)));
% Normalize by trial WITHOUT LED
for i=1:size(newDataY1,1)
    if keepZero==1
        height=mean(newDataY1(i,baselineInds));
    else
        height=mean(newDataY1(i,peakInds))-mean(newDataY1(i,baselineInds));
    end
    if scale==1
        normFactor(i)=1/height;
    else
        normFactor(i)=1;
    end
end
for i=1:size(newDataY1,1)
    if keepZero==1
        newDataY1(i,:)=newDataY1(i,:)*normFactor(i);
        newDataY2(i,:)=newDataY2(i,:)*normFactor(i);
    else
        newDataY1(i,:)=(newDataY1(i,:)-mean(newDataY1(i,baselineInds)))*normFactor(i);
        newDataY2(i,:)=(newDataY2(i,:)-mean(newDataY2(i,baselineInds)))*normFactor(i);
    end
end
% newDataY1=newDataY1+sharedBaseline;
% newDataY2=newDataY2+sharedBaseline;

% Get std. dev. measurements
stds1=std(newDataY1,[],1);
stds2=std(newDataY2,[],1);

% All lines figures
figure();
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2]};
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
% cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
%     [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
%     [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
%     [0.15 0 0.7],[0.7 0.95 0.8],[0.95 0.7 0.8],[0.8 0.65 0.7],[0.35 0.7 0.8],[0.2 0.45 0.2]};
cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};

% color legend
for i=1:length(cs)
    line([i i],[0 1],'Color',cs{i});
end
figure();
for i=1:size(newDataY1,1)
    plot(newXpoints,newDataY1(i,:),'Color',cs{i});
    hold on;
end
figure();
for i=1:size(newDataY2,1)
    plot(newXpoints,newDataY2(i,:),'Color',cs{i});
    hold on;
end

figure();
hax=axes();
hl=plot(newXpoints,mean(newDataY1,1),'Color','k');
addErrBar(newXpoints,mean(newDataY1,1),stds1,'y',hax,hl);
hold on;
hl=plot(newXpoints,mean(newDataY2,1),'Color','r');
addErrBar(newXpoints,mean(newDataY2,1),stds2,'y',hax,hl);

pvals=zeros(length(newXpoints));
for i=1:length(newXpoints)
    [h,p,ci]=ttest(newDataY1(:,i),newDataY2(:,i));
    pvals(i)=p;
end

maxBoth=max(max(mean(newDataY1,1)),max(mean(newDataY2,1)));
minBoth=min(min(mean(newDataY1,1)),min(mean(newDataY2,1)));
halfX=(newXpoints(2)-newXpoints(1))/2;
for i=1:length(newXpoints)
    if pvals(i)<alpha1
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth maxBoth],'Color',[0.5 0.5 0.5],'LineWidth',3);
    end
    above=(maxBoth-minBoth)/15;
    if pvals(i)<alpha2
        line([newXpoints(i)-halfX newXpoints(i)+halfX],[maxBoth+above maxBoth+above],'Color','k','LineWidth',2);
    end
end
end
