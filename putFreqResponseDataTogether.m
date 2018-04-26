function [allTogether_SUmatrix,allTogether_freqTuning,unitFiringRates,allTogether_SUmatrixsingleTrials,finalDiagsOut,MUdiags,unitsOUT,redAndBlue]=putFreqResponseDataTogether(useDir,listing)
%listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');

doAmp=1;
useOnly=0;
useFS=0;
useOnlyLayers=1;
useLayer='L4';
minusOffDiags=0;
firstHarmonic=0;
patchBaseline=0;

finalDiagsOut=[];
MUdiags=[];
% listing=dir(useDir);
% Put PSTH data together
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% allTogether_SUmatrix=zeros(15,15);
% allTogether_SUmatrixsingleTrials=zeros(15,15);
allTogether_freqTuning=zeros(300,15);
allTogether_freqTuningsingleTrials=zeros(300,15);
avEvFR=cell(1,15);
avEvFR_noBaseSub=cell(1,15);
avEvFR_OUT=cell(1,15);
avEvFR_noBaseSub_OUT=cell(1,15);
avEvFR_RS=zeros(15,1);
avEvFR_noBaseSub_RS=zeros(15,1);
avEvFR_FS=zeros(15,1);
avEvFR_noBaseSub_FS=zeros(15,1);
altogetherDiags=[];
altogetherDiagssingleTrials=[];
nUnits=0;
MUdiags=[];
countingMU=1;
countingRS=0;
countingFS=0;
whichGroup=[];
allTogether_SUmatrix=[];
allTogether_SUmatrixsingleTrials=[];
concatEvChs=[];
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
    if ~exist([useDir '\' currFolder '\MUcrosscorr_1.mat'],'file')
        MUcrosscorr=[];
    else
        a=load([useDir '\' currFolder '\MUcrosscorr_1.mat']);
        MUcrosscorr=a.MUcrosscorr;
        if doAmp==1
            MUcrosscorr.pStim{1}=sqrt(MUcrosscorr.pStim{1});
        end
    end
    if ~exist([useDir '\' currFolder '\CSDalignment.mat'],'file')
        CSDalignment=[];
    else
        a=load([useDir '\' currFolder '\CSDalignment.mat']);
        CSDalignment=a.CSDalignment;
    end
    if ~exist([useDir '\' currFolder '\allSU1.mat'],'file')
        allSU=[];
    else
        a=load([useDir '\' currFolder '\allSU1.mat']);
        allSU=a.allSU;
        if doAmp==1
            for i=1:length(allSU)
                allSU(i).pStim{1}=sqrt(allSU(i).pStim{1});
            end
        end
    end
    if ~exist([useDir '\' currFolder '\allSUsingleTrials1.mat'],'file')
        allSUsingleTrials=[];
    else
        a=load([useDir '\' currFolder '\allSUsingleTrials1.mat']);
        allSUsingleTrials=a.allSU;
        if doAmp==1
            for i=1:length(allSUsingleTrials)
                allSUsingleTrials(i).pStim{1}=sqrt(allSUsingleTrials(i).pStim{1});
            end
        end
    end
    if ~exist([useDir '\' currFolder '\criterionUnitsInfo.mat'],'file')
        newassignsinfo=[];
    else
        a=load([useDir '\' currFolder '\criterionUnitsInfo.mat']);
        newassignsinfo=a.newassignsinfo;
    end
    if ~exist([useDir '\' currFolder '\baseUnitsFx.mat'],'file')
        baseUnitsFx=[];
    else
        a=load([useDir '\' currFolder '\baseUnitsFx.mat']);
        baseUnitsFx=a.baseUnitsFx;
    end
    if ~exist([useDir '\' currFolder '\unitsFx.mat'],'file')
        unitsFx=[];
    else
        a=load([useDir '\' currFolder '\unitsFx.mat']);
        unitsFx=a.unitsFx;
    end
    meanBase=zeros(15,1);
    meanEv=zeros(15,1);  
    meanBaseRS=zeros(15,1);
    meanEvRS=zeros(15,1);  
    meanBaseFS=zeros(15,1);
    meanEvFS=zeros(15,1);  
    if ~isempty(MUcrosscorr)
        pStim=MUcrosscorr.pStim{1};
        offDiags=zeros(1,size(pStim,2));
        for j=1:size(pStim,2)
            currRow=j;
            offDiags(j)=nanmean(pStim([1:currRow-1 currRow+1:end],j));
        end
        offDiagMatrix=repmat(offDiags,size(pStim,1),1);
        offDiagMatrix(isnan(offDiagMatrix))=0;
        if patchBaseline==1
            patchVal=sum(sum(pStim(1:3,7:15)))./27;
            pStimMinusOffDiags=pStim-patchVal;
        elseif minusOffDiags==1
            pStimMinusOffDiags=pStim-offDiagMatrix;
        else
            pStimMinusOffDiags=pStim;
        end
        diags=zeros(1,size(pStimMinusOffDiags,1));
        for j=1:length(diags)
            diags(j)=pStimMinusOffDiags(j,j);
        end
%         if doAmp==1
%             MUdiags(countingMU,:)=sqrt(diags);
%         else
            MUdiags(countingMU,:)=diags;
%         end
        countingMU=countingMU+1;
    end
    if ~isempty(baseUnitsFx) && ~isempty(unitsFx)
        % All units
        for j=1:15
            avEvFR{j}=[avEvFR{j} unitsFx.fr{j}-baseUnitsFx.fr{j}];
            avEvFR_noBaseSub{j}=[avEvFR_noBaseSub{j} unitsFx.fr{j}];
        end
        % RS
        [so,so_ind]=sort(newassignsinfo.new_assigns);
        reorder_isFS=newassignsinfo.isFS(so_ind);
        for j=1:15
            curr=baseUnitsFx.fr{j};
            meanBaseRS(j)=nanmean(curr(reorder_isFS==0));
            suBaseRS=curr(reorder_isFS==1);
            curr=unitsFx.fr{j};
            meanEvRS(j)=nanmean(curr(reorder_isFS==0));
            suRS=curr(reorder_isFS==1);
            avEvFR_OUT{j}=[avEvFR_OUT{j} suRS-suBaseRS];
            avEvFR_noBaseSub_OUT{j}=[avEvFR_noBaseSub_OUT{j} suRS];
        end
%         countingRS=countingRS+sum(reorder_isFS==0);
        countingRS=countingRS+1;
        avEvFR_RS=avEvFR_RS+(meanEvRS-meanBaseRS);
        avEvFR_noBaseSub_RS=avEvFR_noBaseSub_RS+meanEvRS;
        % FS
        for j=1:15
            curr=baseUnitsFx.fr{j};
            meanBaseFS(j)=nanmean(curr(reorder_isFS==1));
            curr=unitsFx.fr{j};
            meanEvFS(j)=nanmean(curr(reorder_isFS==1));
        end
        toUseNotNan=~isnan(meanEvFS) & ~isnan(meanBaseFS);
%         countingFS=countingFS+sum(reorder_isFS==1);
        countingFS=countingFS+1;
        avEvFR_FS(toUseNotNan)=avEvFR_FS(toUseNotNan)+(meanEvFS(toUseNotNan)-meanBaseFS(toUseNotNan));
        avEvFR_noBaseSub_FS(toUseNotNan)=avEvFR_noBaseSub_FS(toUseNotNan)+meanEvFS(toUseNotNan);
    end
    newassignsinfo.calibrated_evCh=16-CSDalignment+newassignsinfo.calibrated_evCh;
    concatEvChs=[concatEvChs newassignsinfo.calibrated_evCh];
    if ~isempty(newassignsinfo) && ~isempty(allSU) 
        if useOnly==1
            if useFS==1
                useTheUnits=newassignsinfo.isFS==1;
            else
                useTheUnits=newassignsinfo.isFS==0;
            end
            allSU=allSU(useTheUnits);
            if ~isempty(allSUsingleTrials)
                allSUsingleTrials=allSUsingleTrials(useTheUnits);
            end
        elseif useOnlyLayers==1
            if strcmp(useLayer,'L23')
                useTheUnits=newassignsinfo.calibrated_evCh>=0 & newassignsinfo.calibrated_evCh<=15.84;
            elseif strcmp(useLayer,'L4')
                useTheUnits=newassignsinfo.calibrated_evCh>15.84 & newassignsinfo.calibrated_evCh<=17.4432;
            elseif strcmp(useLayer,'L5a')
                useTheUnits=newassignsinfo.calibrated_evCh>17.4432 & newassignsinfo.calibrated_evCh<=18.24;
            elseif strcmp(useLayer,'L5b')
                useTheUnits=newassignsinfo.calibrated_evCh>18.24 & newassignsinfo.calibrated_evCh<=19.8348;
            elseif strcmp(useLayer,'L6')
                useTheUnits=newassignsinfo.calibrated_evCh>19.8348;
            end
            allSU=allSU(useTheUnits);
            if ~isempty(allSUsingleTrials)
                allSUsingleTrials=allSUsingleTrials(useTheUnits);
            end
        end                
        altogetherDiagssingleTrials=[];
        for j=1:length(allSU)
            nUnits=nUnits+1;
            currSU=allSU(j);
            pStim=currSU.pStim{1};
            if ~isempty(allSUsingleTrials)
                currSUsingleTrials=allSUsingleTrials(j);
                pStimsingleTrials=currSUsingleTrials.pStim{1};
            end
%             if doAmp==1
%                 pStim=sqrt(pStim);
%                 if ~isempty(allSUsingleTrials)
%                     pStimsingleTrials=sqrt(pStimsingleTrials);
%                 end
%             end
%             if any(any(isnan(pStim)))==1
%                 pStim(isnan(pStim))=0;
%             end
%             if any(any(isnan(pStimsingleTrials)))==1
%                 pStimsingleTrials(isnan(pStimsingleTrials))=0;
%             end
%             if doAmp==1
%                 allTogether_SUmatrix(nUnits,:,:)=sqrt(pStim);
%             else
                allTogether_SUmatrix(nUnits,:,:)=pStim;
%             end
%             allTogether_SUmatrix{nUnits}=pStim;
            offDiags=zeros(1,size(pStim,2));
            if ~isempty(allSUsingleTrials)
                offDiagssingleTrials=zeros(1,size(pStimsingleTrials,2));
            end
            for k=1:size(pStim,2)
                currRow=k;
                nonspec=pStim([1:currRow-1 currRow+1:end],k);
                offDiags(k)=mean(nonspec(~isnan(nonspec)));
            end
            if ~isempty(allSUsingleTrials)
                for k=1:size(pStimsingleTrials,2)
                    currRow=k;
                    nonspec=pStimsingleTrials([1:currRow-1 currRow+1:end],k);
                    offDiagssingleTrials(k)=mean(nonspec(~isnan(nonspec)));
                end
            end
            offDiagMatrix=repmat(offDiags,size(pStim,1),1);
            offDiagMatrix(isnan(offDiagMatrix))=0;
            if patchBaseline==1
                patchVal=sum(sum(pStim(1:3,7:15)))./27;
                pStimMinusOffDiags=pStim-patchVal;
            elseif minusOffDiags==1
                pStimMinusOffDiags=pStim-offDiagMatrix;
            else
                pStimMinusOffDiags=pStim;
            end
            diags=zeros(1,size(pStimMinusOffDiags,1)); 
            if ~isempty(allSUsingleTrials)
                offDiagMatrixsingleTrials=repmat(offDiagssingleTrials,size(pStimsingleTrials,1),1);
                offDiagMatrixsingleTrials(isnan(offDiagMatrixsingleTrials))=0;
                pStimMinusOffDiagssingleTrials=pStimsingleTrials-offDiagMatrixsingleTrials;
                allTogether_SUmatrixsingleTrials(nUnits,:,:)=pStimMinusOffDiagssingleTrials;
                diagssingleTrials=zeros(1,size(pStimMinusOffDiagssingleTrials,1));
            end
            if firstHarmonic==1 
                diags(1)=pStimMinusOffDiags(1,2);
                diags(2)=pStimMinusOffDiags(2,3);
                diags(3)=pStimMinusOffDiags(3,5);
                diags(4)=pStimMinusOffDiags(4,7);
                diags(5)=pStimMinusOffDiags(5,9);
                diags(6)=pStimMinusOffDiags(6,11);
                diags(7)=nan;
                diags(8)=nan;
                diags(9)=nan;
                diags(10)=nan;
                diags(11)=pStimMinusOffDiags(11,13);
                diags(12)=pStimMinusOffDiags(12,15);
                diags(13)=nan;
                diags(14)=nan;
                diags(15)=nan;
            else
                for k=1:length(diags)
                    diags(k)=pStimMinusOffDiags(k,k);
                end
            end
            if ~isempty(allSUsingleTrials)
                for k=1:length(diagssingleTrials)
                    diagssingleTrials(k)=pStimMinusOffDiagssingleTrials(k,k);
                end
            end
            currAss=currSU.assigns{1};
            calEvCh=newassignsinfo.calibrated_evCh(newassignsinfo.new_assigns==currAss);
%             if any(isnan(diags))==1
%                 disp('problem');
%             end
%             if doAmp==1
%                 diags=sqrt(diags);
%                 if ~isempty(allSUsingleTrials)
%                     diagssingleTrials=sqrt(diagssingleTrials);
%                 end
%             end
            altogetherDiags(nUnits,:)=diags;
            whichGroup(nUnits)=floor(calEvCh/1.5);
            allTogether_freqTuning(16-CSDalignment+floor(calEvCh/3),:)=allTogether_freqTuning(16-CSDalignment+floor(calEvCh/3),:)+diags;
            if ~isempty(allSUsingleTrials)
                altogetherDiagssingleTrials(nUnits,:)=diagssingleTrials;
                allTogether_freqTuningsingleTrials(16-CSDalignment+floor(calEvCh/3),:)=allTogether_freqTuningsingleTrials(16-CSDalignment+floor(calEvCh/3),:)+diagssingleTrials;
            end  
        end
    end
end
for i=1:15
    avEvFRmean(i)=nanmean(avEvFR{i});
    avEvFR_noBaseSubmean(i)=nanmean(avEvFR_noBaseSub{i});
end
disp('avEvFR');
disp(avEvFRmean);
disp('avEvFR_noBaseSub');
disp(avEvFR_noBaseSubmean);
unitFiringRates.avEvFR=avEvFR;
unitFiringRates.avEvFR_noBaseSub=avEvFR_noBaseSub;
unitsOUT.avEvFR=avEvFR_OUT;
unitsOUT.avEvFR_noBaseSub=avEvFR_noBaseSub_OUT;

avEvFR_noBaseSub_RS=avEvFR_noBaseSub_RS/countingRS;
avEvFR_RS=avEvFR_RS/countingRS;
avEvFR_noBaseSub_FS=avEvFR_noBaseSub_FS/countingFS;
avEvFR_FS=avEvFR_FS/countingFS;
figure(); 
plot(freqs,avEvFR_noBaseSub_RS,'Color','r');
hold on; 
plot(freqs,avEvFR_noBaseSub_FS,'Color','b');
redAndBlue.RS=avEvFR_noBaseSub_RS;
redAndBlue.FS=avEvFR_noBaseSub_FS;
title('RS evoked red, FS evoked blue');

if ~isempty(MUdiags)
    figure(); 
    semilogx(freqs,nanmean(MUdiags,1));
    title('Amp of MU freq response');
end

if ~isempty(altogetherDiags)
    figure(); 
    plot(freqs,altogetherDiags','Color',[0.7 0.7 0.7]);
    hold on;
    plot(freqs,nanmean(altogetherDiags,1),'Color','k');
    addErrBar(freqs,nanmean(altogetherDiags,1),nanstd(altogetherDiags,[],1)./sqrt(size(altogetherDiags,1)),'y',[]);
    title('Average Freq. Response for Average Unit PSTH');
%     finalDiagsOut=mean(altogetherDiags,1);
    finalDiagsOut=altogetherDiags;

    norm_altogetherDiags=altogetherDiags;
    thePeaks=nanmax(altogetherDiags,[],2);
    scale=repmat(thePeaks,1,size(altogetherDiags,2));
    norm_altogetherDiags=norm_altogetherDiags./scale;
%     finalDiagsOut=norm_altogetherDiags;

    
    figure();
    plot(freqs,norm_altogetherDiags','Color',[0.7 0.7 0.7]);
    hold on;
    plot(freqs,nanmean(norm_altogetherDiags,1),'Color','k');
    addErrBar(freqs,nanmean(norm_altogetherDiags,1),nanstd(norm_altogetherDiags,[],1)./sqrt(size(norm_altogetherDiags,1)),'y',[]);
    title('Average Freq. Response for Average Unit PSTH -- NORMED');
    
    figure();
    w=unique(whichGroup);
    for i=1:length(w)
        currGroup=w(i);
        plot(freqs,nanmean(norm_altogetherDiags(whichGroup==currGroup,:),1),'Color',cs{i});
        addErrBar(freqs,nanmean(norm_altogetherDiags(whichGroup==currGroup,:),1),nanstd(norm_altogetherDiags(whichGroup==currGroup,:),[],1)./sqrt(size(norm_altogetherDiags(whichGroup==currGroup,:),1)),'y',[],[],cs{i});
        hold all;
    end
end

if ~isempty(altogetherDiagssingleTrials)
    figure();
    plot(freqs,altogetherDiagssingleTrials','Color',[0.7 0.7 0.7]);
    hold on;
    plot(freqs,nanmean(altogetherDiagssingleTrials,1),'Color','k');
    addErrBar(freqs,nanmean(altogetherDiagssingleTrials,1),nanstd(altogetherDiagssingleTrials,[],1)./sqrt(size(altogetherDiagssingleTrials,1)),'y',[]);
    title('Average Freq. Response for Unit Single Trials');
%     finalDiagsOut=nanmean(altogetherDiagssingleTrials,1);

    norm_altogetherDiagssingleTrials=altogetherDiagssingleTrials;
    thePeaks=nanmax(altogetherDiagssingleTrials,[],2);
    scale=repmat(thePeaks,1,size(altogetherDiagssingleTrials,2));
    norm_altogetherDiagssingleTrials=norm_altogetherDiagssingleTrials./scale;

    figure();
    plot(freqs,norm_altogetherDiagssingleTrials','Color',[0.7 0.7 0.7]);
    hold on;
    plot(freqs,nanmean(norm_altogetherDiagssingleTrials,1),'Color','k');
    addErrBar(freqs,nanmean(norm_altogetherDiagssingleTrials,1),nanstd(norm_altogetherDiagssingleTrials,[],1)./sqrt(size(norm_altogetherDiagssingleTrials,1)),'y',[]);
    title('Average Freq. Response for Unit Single Trials -- NORMED');
end
    
if ~isempty(allTogether_SUmatrix)
    figure(); 
    imagesc(reshape(nanmean(allTogether_SUmatrix,1),[15 15]));
    title('Average PSTH SU Matrix');
end
if ~isempty(allTogether_SUmatrixsingleTrials)
    figure();
    imagesc(reshape(nanmean(allTogether_SUmatrixsingleTrials,1),[15 15]));
    title('Single Trial SU Matrix');
end

% if ~isempty(allTogether_freqTuning)
%     figure();
%     imagesc(allTogether_freqTuning);
%     title('Average PSTH Tuning');
%     normTuning=allTogether_freqTuning;
%     for i=1:size(normTuning,1)
%         mi=min(normTuning(i,:));
%         normTuning(i,:)=normTuning(i,:)-mi;
%         ma=max(allTogether_freqTuning(i,:));
%         if ma>0
%             normTuning(i,:)=normTuning(i,:)/ma;
%         end
%     end
%     figure();
%     imagesc(normTuning);
%     title('Average PSTH Tuning -- NORMED');
%     
%     normTuning=allTogether_freqTuningsingleTrials;
%     for i=1:size(normTuning,1)
%         mi=min(normTuning(i,:));
%         normTuning(i,:)=normTuning(i,:)-mi;
%         ma=max(allTogether_freqTuningsingleTrials(i,:));
%         if ma>0
%             normTuning(i,:)=normTuning(i,:)/ma;
%         end
%     end
%     figure();
%     imagesc(normTuning);
%     title('Single Trial Tuning -- NORMED');
% end
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
