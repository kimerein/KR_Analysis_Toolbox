function putFreqLayersTogether(useDir,listing)
%listing=dir('W:\Analysis Computer\Thalamus Silencing Across Mice\ANESTH\Evoked dLGN');

doAmp=1;
useOnly=0;
useFS=1;
useOnlyLayers=1;
useLayer='L5a';

finalDiagsOut=[];
MUdiags=[];
% listing=dir(useDir);
% Put PSTH data together
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% allTogether_SUmatrix=zeros(15,15);
% allTogether_SUmatrixsingleTrials=zeros(15,15);
allTogether_freqTuning=zeros(300,15);
allTogether_freqTuningsingleTrials=zeros(300,15);
avEvFR=zeros(15,1);
avEvFR_noBaseSub=zeros(15,1);
avEvFR_RS=zeros(15,1);
avEvFR_noBaseSub_RS=zeros(15,1);
avEvFR_FS=zeros(15,1);
avEvFR_noBaseSub_FS=zeros(15,1);
altogetherDiags=[];
altogetherDiagssingleTrials=[];
nUnits=0;
MUdiags=[];
MUlayers=zeros(300,15);
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
    if ~exist([useDir '\' currFolder '\MUcrosscorr_1.mat'],'file')
        MUcrosscorr=[];
    else
        a=load([useDir '\' currFolder '\MUcrosscorr_1.mat']);
        MUcrosscorr=a.MUcrosscorr;
    end
    if ~exist([useDir '\' currFolder '\CSDalignment.mat'],'file')
        CSDalignment=[];
    else
        a=load([useDir '\' currFolder '\CSDalignment.mat']);
        CSDalignment=a.CSDalignment;
    end
    if ~isempty(MUcrosscorr)
        pStim=MUcrosscorr.pStim{1};
        if doAmp==1
            pStim=sqrt(pStim);
        end
        offDiags=zeros(1,size(pStim,2));
        for j=1:size(pStim,2)
            currRow=j;
            offDiags(j)=nanmean(pStim([1:currRow-1 currRow+1:end],j));
        end
        offDiagMatrix=repmat(offDiags,size(pStim,1),1);
        offDiagMatrix(isnan(offDiagMatrix))=0;
        pStimMinusOffDiags=pStim-offDiagMatrix;
        diags=zeros(1,size(pStimMinusOffDiags,1));
        for j=1:length(diags)
            diags(j)=pStimMinusOffDiags(j,j);
        end
        MUdiags(countingMU,:)=diags;
        countingMU=countingMU+1;
    end
    MUlayers(150+16-CSDalignment:150+16-CSDalignment-1+size(MUdiags,1)       
        altogetherDiagssingleTrials=[];
        for j=1:length(allSU)
            nUnits=nUnits+1;
            currSU=allSU(j);
            pStim=currSU.pStim{1};
            if ~isempty(allSUsingleTrials)
                currSUsingleTrials=allSUsingleTrials(j);
                pStimsingleTrials=currSUsingleTrials.pStim{1};
            end
            if doAmp==1
                pStim=sqrt(pStim);
                if ~isempty(allSUsingleTrials)
                    pStimsingleTrials=sqrt(pStimsingleTrials);
                end
            end
%             if any(any(isnan(pStim)))==1
%                 pStim(isnan(pStim))=0;
%             end
%             if any(any(isnan(pStimsingleTrials)))==1
%                 pStimsingleTrials(isnan(pStimsingleTrials))=0;
%             end
            allTogether_SUmatrix(nUnits,:,:)=pStim;
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
            pStimMinusOffDiags=pStim-offDiagMatrix;
            diags=zeros(1,size(pStimMinusOffDiags,1)); 
            if ~isempty(allSUsingleTrials)
                offDiagMatrixsingleTrials=repmat(offDiagssingleTrials,size(pStimsingleTrials,1),1);
                offDiagMatrixsingleTrials(isnan(offDiagMatrixsingleTrials))=0;
                pStimMinusOffDiagssingleTrials=pStimsingleTrials-offDiagMatrixsingleTrials;
                allTogether_SUmatrixsingleTrials(nUnits,:,:)=pStimMinusOffDiagssingleTrials;
                diagssingleTrials=zeros(1,size(pStimMinusOffDiagssingleTrials,1));
            end
            for k=1:length(diags)
                diags(k)=pStimMinusOffDiags(k,k);
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
avEvFR=avEvFR/nUnits;
avEvFR_noBaseSub=avEvFR_noBaseSub/nUnits;
disp('avEvFR');
disp(avEvFR);
disp('avEvFR_noBaseSub');
disp(avEvFR_noBaseSub);

avEvFR_noBaseSub_RS=avEvFR_noBaseSub_RS/countingRS;
avEvFR_RS=avEvFR_RS/countingRS;
avEvFR_noBaseSub_FS=avEvFR_noBaseSub_FS/countingFS;
avEvFR_FS=avEvFR_FS/countingFS;
figure(); 
plot(freqs,avEvFR_noBaseSub_RS,'Color','r');
hold on; 
plot(freqs,avEvFR_noBaseSub_FS,'Color','b');
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
