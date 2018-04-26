function [allTogether_SUmatrix,allTogether_freqTuning,avEvFR_noBaseSub,allTogether_SUmatrixsingleTrials,finalDiagsOut]=putFreqResponseDataTogether(useDir,listing)
%listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');

doAmp=1;
useOnly=1;
useFS=0;

% listing=dir(useDir);
% Put PSTH data together
gotData=0;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
allTogether_SUmatrix=zeros(15,15);
allTogether_SUmatrixsingleTrials=zeros(15,15);
allTogether_freqTuning=zeros(300,15);
allTogether_freqTuningsingleTrials=zeros(300,15);
avEvFR=zeros(15,1);
avEvFR_noBaseSub=zeros(15,1);
altogetherDiags=[];
altogetherDiagssingleTrials=[];
nUnits=0;
whichGroup=[];
cs={[0 0 0],[0 0.8 0.8],[0.8 0 0.8],[0.8 0.8 0],[0.2 0.2 0.8],[0.2 0 0.2],...
    [0.5 0 0],[0.5 0.8 0.8],[0.8 0.5 0.8],[0.8 0.8 0.5],[0.2 0.5 0.8],[0.2 0.5 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2],...
    [0 0 0.7],[0.7 0.8 0.8],[0.8 0.7 0.8],[0.8 0.8 0.7],[0.2 0.7 0.8],[0.2 0.7 0.2]};
for i=1:length(listing)
    currFolder=listing(i).name;
%     if ~exist([useDir '\' currFolder '\exptData.mat'],'file')
%         continue
%     end
%     a=load([useDir '\' currFolder '\exptData.mat']);
    if ~exist([useDir '\' currFolder '\CSDalignment.mat'],'file')
        continue
    end
    gotData=gotData+1;
    a=load([useDir '\' currFolder '\CSDalignment.mat']);
    CSDalignment=a.CSDalignment;
    a=load([useDir '\' currFolder '\allSU1.mat']);
    allSU=a.allSU;
    a=load([useDir '\' currFolder '\allSUsingleTrials1.mat']);
    allSUsingleTrials=a.allSU;
    a=load([useDir '\' currFolder '\criterionUnitsInfo.mat']);
    newassignsinfo=a.newassignsinfo;
    a=load([useDir '\' currFolder '\baseUnitsFx.mat']);
    baseUnitsFx=a.baseUnitsFx;
    a=load([useDir '\' currFolder '\unitsFx.mat']);
    unitsFx=a.unitsFx;
    meanBase=zeros(15,1);
    meanEv=zeros(15,1);    
    for j=1:15
        meanBase(j)=mean(baseUnitsFx.fr{j});
        meanEv(j)=mean(unitsFx.fr{j});
    end
    avEvFR=avEvFR+(meanEv-meanBase);
    avEvFR_noBaseSub=avEvFR_noBaseSub+meanEv;
    if useOnly==1
        if useFS==1
            useTheUnits=newassignsinfo.isFS==1;
        else
            useTheUnits=newassignsinfo.isFS==0;
        end
        allSU=allSU(useTheUnits);
        allSUsingleTrials=allSUsingleTrials(useTheUnits);
    end
    for j=1:length(allSU)
        nUnits=nUnits+1;
        currSU=allSU(j);
        currSUsingleTrials=allSUsingleTrials(j);
        pStim=currSU.pStim{1};
        pStimsingleTrials=currSUsingleTrials.pStim{1};
        if doAmp==1
            pStim=sqrt(pStim);
            pStimsingleTrials=sqrt(pStimsingleTrials);
        end
        if any(any(isnan(pStim)))==1
            pStim(isnan(pStim))=0;
        end
        if any(any(isnan(pStimsingleTrials)))==1
            pStimsingleTrials(isnan(pStimsingleTrials))=0;
        end
        allTogether_SUmatrix=allTogether_SUmatrix+pStim;
        offDiags=zeros(1,size(pStim,2));
        offDiagssingleTrials=zeros(1,size(pStimsingleTrials,2));
        for k=1:size(pStim,2)
            currRow=k;
            offDiags(k)=mean(pStim([1:currRow-1 currRow+1:end],k));
        end
        for k=1:size(pStimsingleTrials,2)
            currRow=k;
            offDiagssingleTrials(k)=mean(pStimsingleTrials([1:currRow-1 currRow+1:end],k));
        end
        offDiagMatrix=repmat(offDiags,size(pStim,1),1);
        offDiagMatrixsingleTrials=repmat(offDiagssingleTrials,size(pStimsingleTrials,1),1);
%         offDiagMatrixsingleTrials=repmat(offDiagssingleTrials,1,size(pStimsingleTrials,2));
%         offDiagMatrixsingleTrials=zeros(size(offDiagMatrix));
        pStimMinusOffDiags=pStim-offDiagMatrix;
        pStimMinusOffDiagssingleTrials=pStimsingleTrials-offDiagMatrixsingleTrials;
        allTogether_SUmatrixsingleTrials=allTogether_SUmatrixsingleTrials+pStimMinusOffDiagssingleTrials;
        diags=zeros(1,size(pStimMinusOffDiags,1));
        diagssingleTrials=zeros(1,size(pStimMinusOffDiagssingleTrials,1));
        for k=1:length(diags)
            diags(k)=pStimMinusOffDiags(k,k);
        end
        for k=1:length(diagssingleTrials)
            diagssingleTrials(k)=pStimMinusOffDiagssingleTrials(k,k);
        end
        currAss=currSU.assigns{1};
        calEvCh=newassignsinfo.calibrated_evCh(newassignsinfo.new_assigns==currAss);
        if any(isnan(diags))==1
            disp('problem');
        end
        altogetherDiags(nUnits,:)=diags;
        altogetherDiagssingleTrials(nUnits,:)=diagssingleTrials;
        whichGroup(nUnits)=floor(calEvCh/1.5);
        allTogether_freqTuning(16-CSDalignment+floor(calEvCh/3),:)=allTogether_freqTuning(16-CSDalignment+floor(calEvCh/3),:)+diags;
        allTogether_freqTuningsingleTrials(16-CSDalignment+floor(calEvCh/3),:)=allTogether_freqTuningsingleTrials(16-CSDalignment+floor(calEvCh/3),:)+diagssingleTrials;
    end
end
allTogether_SUmatrix=allTogether_SUmatrix/nUnits;
allTogether_SUmatrixsingleTrials=allTogether_SUmatrixsingleTrials/nUnits;
avEvFR=avEvFR/nUnits;
avEvFR_noBaseSub=avEvFR_noBaseSub/nUnits;
disp('avEvFR');
disp(avEvFR);
disp('avEvFR_noBaseSub');
disp(avEvFR_noBaseSub);

figure(); 
plot(freqs,altogetherDiags','Color',[0.7 0.7 0.7]);
hold on;
plot(freqs,mean(altogetherDiags,1),'Color','k');
addErrBar(freqs,mean(altogetherDiags,1),std(altogetherDiags,[],1)./sqrt(size(altogetherDiags,1)),'y',[]);
title('Average Freq. Response for Average Unit PSTH');
% finalDiagsOut=mean(altogetherDiags,1);

norm_altogetherDiags=altogetherDiags;
thePeaks=max(altogetherDiags,[],2);
scale=repmat(thePeaks,1,size(altogetherDiags,2));
norm_altogetherDiags=norm_altogetherDiags./scale;

figure(); 
plot(freqs,norm_altogetherDiags','Color',[0.7 0.7 0.7]);
hold on;
plot(freqs,mean(norm_altogetherDiags,1),'Color','k');
addErrBar(freqs,mean(norm_altogetherDiags,1),std(norm_altogetherDiags,[],1)./sqrt(size(norm_altogetherDiags,1)),'y',[]);
title('Average Freq. Response for Average Unit PSTH -- NORMED');

figure(); 
plot(freqs,altogetherDiagssingleTrials','Color',[0.7 0.7 0.7]);
hold on;
plot(freqs,mean(altogetherDiagssingleTrials,1),'Color','k');
addErrBar(freqs,mean(altogetherDiagssingleTrials,1),std(altogetherDiagssingleTrials,[],1)./sqrt(size(altogetherDiagssingleTrials,1)),'y',[]);
title('Average Freq. Response for Unit Single Trials');
finalDiagsOut=mean(altogetherDiagssingleTrials,1);

norm_altogetherDiagssingleTrials=altogetherDiagssingleTrials;
thePeaks=max(altogetherDiagssingleTrials,[],2);
scale=repmat(thePeaks,1,size(altogetherDiagssingleTrials,2));
norm_altogetherDiagssingleTrials=norm_altogetherDiagssingleTrials./scale;

figure(); 
plot(freqs,norm_altogetherDiagssingleTrials','Color',[0.7 0.7 0.7]);
hold on;
plot(freqs,mean(norm_altogetherDiagssingleTrials,1),'Color','k');
addErrBar(freqs,mean(norm_altogetherDiagssingleTrials,1),std(norm_altogetherDiagssingleTrials,[],1)./sqrt(size(norm_altogetherDiagssingleTrials,1)),'y',[]);
title('Average Freq. Response for Unit Single Trials -- NORMED');

figure(); 
w=unique(whichGroup);
for i=1:length(w)
    currGroup=w(i);
    plot(freqs,mean(norm_altogetherDiags(whichGroup==currGroup,:),1),'Color',cs{i});
    addErrBar(freqs,mean(norm_altogetherDiags(whichGroup==currGroup,:),1),std(norm_altogetherDiags(whichGroup==currGroup,:),[],1)./sqrt(size(norm_altogetherDiags(whichGroup==currGroup,:),1)),'y',[],[],cs{i});
    hold all;
end

figure(); 
imagesc(allTogether_SUmatrix);
title('Average PSTH SU Matrix');
figure(); 
imagesc(allTogether_SUmatrixsingleTrials);
title('Single Trial SU Matrix');
figure(); 
imagesc(allTogether_freqTuning);
title('Average PSTH Tuning');

normTuning=allTogether_freqTuning;
for i=1:size(normTuning,1)
    mi=min(normTuning(i,:));
    normTuning(i,:)=normTuning(i,:)-mi;
    ma=max(allTogether_freqTuning(i,:));
    if ma>0
        normTuning(i,:)=normTuning(i,:)/ma;
    end
end
figure(); 
imagesc(normTuning);
title('Average PSTH Tuning -- NORMED');

normTuning=allTogether_freqTuningsingleTrials;
for i=1:size(normTuning,1)
    mi=min(normTuning(i,:));
    normTuning(i,:)=normTuning(i,:)-mi;
    ma=max(allTogether_freqTuningsingleTrials(i,:));
    if ma>0
        normTuning(i,:)=normTuning(i,:)/ma;
    end
end
figure(); 
imagesc(normTuning);
title('Single Trial Tuning -- NORMED');

% figure();
% scatter(unitEvChs,unitTaus);
% 
% normLayers.y1=togetherLayers.y1;
% normLayers.y2=togetherLayers.y2;
% for i=1:size(normLayers.y1,1)
%     mi=min(normLayers.y1(i,:));
%     normLayers.y1(i,:)=normLayers.y1(i,:)-mi;
%     ma=max(normLayers.y1(i,:));
%     if ma>0
%         normLayers.y1(i,:)=normLayers.y1(i,:)/ma;
%     end
% end
% for i=1:size(normLayers.y2,1)
%     mi=min(normLayers.y2(i,:));
%     normLayers.y2(i,:)=normLayers.y2(i,:)-mi;
%     ma=max(normLayers.y2(i,:));
%     if ma>0
%         normLayers.y2(i,:)=normLayers.y2(i,:)/ma;
%     end
% end
% figure(); 
% imagesc(normLayers.y1);
% figure(); 
% imagesc(normLayers.y2);

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
