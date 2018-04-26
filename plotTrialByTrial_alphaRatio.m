function [sumSpectra,sortedSpectra,bestt,usefinterp]=plotTrialByTrial_alphaRatio(psth,autocorr_out,currUnitMean,psth2,autocorr_out2,evBase,ffrLEDadjusted)

plotNonnormSpecs=1;
normEachTimePoint=0;
regressOffMeans=0;
doCluster=1;
fromPSTH=1;
sortedSpectra=[];

tryS=12;

if length(tryS)==1
    e=evBase(:,tryS);
    [~,evInd]=sort(e);
    [~,ffrInd]=sort(ffrLEDadjusted(:,tryS));
    combInd=zeros(1,length(evInd));
    for i=1:length(evInd)
        combInd(i)=nanmean([find(i==evInd) find(i==ffrInd)]);
    end
    getn=4;
    [~,sortInd]=sort(combInd);
    sortInd=ffrInd';
%     unitsLike=sortInd(end-(getn-1):end);
%     unitsLike=evInd(end-(getn-1):end);
    % unitsLike=ffrInd(end-(getn-1):end);
%     unitsLike=1:13;
    disp(sortInd);
end

% finterp=3.125:0.1:18.75;
% finterp=3:0.1:20-0.1;
finterp=3:0.1:30-0.1;
% finterp=20:0.1:40-0.1;
lagsWindow=[0 2000];
usel=[0];
for i=1:length(psth.psths)
    for z=1:length(tryS)
        uses=tryS(z);
        [f,t,S]=scriptForPlottingAutocorr_sub(psth,autocorr_out,usel,uses,i,lagsWindow,fromPSTH);
        bestt=t;
        for j=1:size(S,1)
            currS(j,:)=interp1(f,S(j,:),finterp);
        end
        autocorr_spectraFirst{z}=currS;
        usefinterp=finterp(~any(isnan(currS),1));
    end
    allS{i}=autocorr_spectraFirst;
end

if ~isempty(psth2)
    lagsWindow=[0 3500];
    usel=[0];
    for i=1:length(psth2.psths)
        for z=1:length(tryS)
            uses=tryS(z);
            [f,t,S]=scriptForPlottingAutocorr_sub(psth2,autocorr_out2,usel,uses,i,lagsWindow,fromPSTH);
            for j=1:size(S,1)
                currS2(j,:)=interp1(f,S(j,:),finterp);
            end
            autocorr_spectraSecond{z}=currS2;
        end
        allSsecond{i}=autocorr_spectraSecond;
    end
end

if regressOffMeans==1
    if isempty(psth2)
        meanForUnits=cell(1,length(psth.psths));
        for i=1:length(psth.psths)
            acrossStimSpectra=[];
            for z=1:length(tryS)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                acrossStimSpectra=[acrossStimSpectra; tempSpectra];
            end
            meanForUnits{i}=nanmean(acrossStimSpectra,1);
        end
        for i=1:length(psth.psths)
            for z=1:length(tryS)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                tempSpectra=tempSpectra-repmat(meanForUnits{i},size(tempSpectra,1),1);
                autocorr_spectraFirst{z}=tempSpectra;
            end
            allS{i}=autocorr_spectraFirst;
        end
    else
        meanForUnits=cell(1,length(psth.psths));
        for i=1:length(psth.psths)
            acrossStimSpectra=[];
            for z=1:length(tryS)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                tempUnit2=allSsecond{i};
                tempSpectra2=tempUnit2{z};
                acrossStimSpectra=[acrossStimSpectra; tempSpectra; tempSpectra2];
            end
            meanForUnits{i}=nanmean(acrossStimSpectra,1);
        end
        for i=1:length(psth.psths)
            for z=1:length(tryS)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                tempUnit2=allSsecond{i};
                tempSpectra2=tempUnit2{z};
                tempSpectra=tempSpectra-repmat(meanForUnits{i},size(tempSpectra,1),1);
                tempSpectra2=tempSpectra2-repmat(meanForUnits{i},size(tempSpectra2,1),1);
                if normEachTimePoint==1
                    for j=1:size(tempSpectra,1)
                        tempSpectra(j,:)=tempSpectra(j,:)-min(tempSpectra(j,:));
                        tempSpectra(j,:)=tempSpectra(j,:)./max(tempSpectra(j,:));
                        tempSpectra2(j,:)=tempSpectra2(j,:)-min(tempSpectra2(j,:));
                        tempSpectra2(j,:)=tempSpectra2(j,:)./max(tempSpectra2(j,:));
                    end
                end
                autocorr_spectraFirst{z}=tempSpectra;
                autocorr_spectraSecond{z}=tempSpectra2;
            end
            allS{i}=autocorr_spectraFirst;
            allSsecond{i}=autocorr_spectraSecond;
        end
    end
end

if plotNonnormSpecs==1
    figure();
    if isempty(psth2)
        for z=1:length(tryS)
            concatSpectra=[];
            sepSpectra=cell(1,length(psth.psths));
            for i=1:length(psth.psths)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                if normEachTimePoint==1
                    for j=1:size(tempSpectra,1)
                        tempSpectra(j,:)=tempSpectra(j,:)-min(tempSpectra(j,:));
                        tempSpectra(j,:)=tempSpectra(j,:)./max(tempSpectra(j,:));
                    end
                end
                concatSpectra=[concatSpectra; tempSpectra'];
                sepSpectra{i}=tempSpectra';
                if i==1
                    sumSpectra=zeros(size(tempSpectra'));
                    sumSpectra_bottom=zeros(size(tempSpectra'));
                    sumSpectra_top=zeros(size(tempSpectra'));
                end
                sumSpectra=sumSpectra+tempSpectra';
                if ismember(i,sortInd(1:getn))
                    sumSpectra_bottom=sumSpectra_bottom+tempSpectra';
                elseif ismember(i,sortInd(end-getn+1:end))
                    sumSpectra_top=sumSpectra_top+tempSpectra';
                end
            end
            subplot(1,length(tryS),z);
            sumSpectra=sumSpectra./length(psth.psths);
            sumSpectra_bottom=sumSpectra_bottom./getn;
            sumSpectra_top=sumSpectra_top./getn;
            heatmap([concatSpectra; sumSpectra]);
            imagesc(bestt,f,[concatSpectra; sumSpectra]);
%             set(gca,'XTick',[1:5:size(concatSpectra,2)]);
%             set(gca,'XTick',[1 floor(size(concatSpectra,2)/2) size(concatSpectra,2)]);
%             xticklab=cell(1,length(1:10:length(t)));
%             k=1;
%             for j=1:10:length(t)
%                 xticklab{k}=num2str(t(j)); % in s
%                 k=k+1;
%             end
%             set(gca,'XTickLabel',{0, 1.5, 3});
            if ~isempty(sortInd)
                figure();
                imagesc(bestt,f,[sumSpectra_bottom; sumSpectra_top]);
%                 set(gca,'XTick',[1:5:size(concatSpectra,2)]);
%                 set(gca,'XTickLabel',{0, 1.5, 3});
                title('Bottom Above, Top Below');
                sortedSpectra=[sumSpectra_bottom; sumSpectra_top];
%                 figure();
%                 heatmap(sumSpectra_top);
%                 set(gca,'XTick',[1 floor(size(concatSpectra,2)/2) size(concatSpectra,2)]);
%                 set(gca,'XTickLabel',{0, 1.5, 3});
%                 title('Top');
                catSpectra=[];
                for i=1:length(sepSpectra)
                    currSpec=sortInd(i);
                    catSpectra=[catSpectra; sepSpectra{currSpec}];
                end
                figure();
                imagesc(bestt,f,catSpectra);
%                 set(gca,'XTick',[1:5:size(concatSpectra,2)]);
%                 set(gca,'XTickLabel',{0,0.25,0.5,0.75,1,1.25,1.5,1.75,2});
%                 set(gca,'XTick',[1 floor(size(concatSpectra,2)/2) size(concatSpectra,2)]);
%                 set(gca,'XTickLabel',{0, 1.5, 3});
            end
            
        end
    else
        for z=1:length(tryS)
            concatSpectra=[];
            for i=1:length(psth.psths)
                tempUnit=allS{i};
                tempSpectra=tempUnit{z};
                tempUnit2=allSsecond{i};
                tempSpectra2=tempUnit2{z};
%                 if normEachTimePoint==1
%                     for j=1:size(tempSpectra,1)
%                         tempSpectra(j,:)=tempSpectra(j,:)-min(tempSpectra(j,:));
%                         tempSpectra(j,:)=tempSpectra(j,:)./max(tempSpectra(j,:));
%                         tempSpectra2(j,:)=tempSpectra2(j,:)-min(tempSpectra2(j,:));
%                         tempSpectra2(j,:)=tempSpectra2(j,:)./max(tempSpectra2(j,:));
%                     end
%                 end
                concatSpectra=[concatSpectra; [tempSpectra' tempSpectra2']];
                if i==1
                    sumSpectra=zeros(size([tempSpectra' tempSpectra2']));
                end
                sumSpectra=sumSpectra+[tempSpectra' tempSpectra2'];
            end
            subplot(1,length(tryS),z);
            sumSpectra=sumSpectra./length(psth.psths);
            heatmap([concatSpectra; sumSpectra]);
            set(gca,'XTick',1:10*2:size(concatSpectra,2));
            xticklab=cell(1,length(1:10:length(t)));
            k=1;
            for j=1:10:length(t)
                xticklab{k}=num2str(t(j)); % in s
                k=k+1;
            end
            set(gca,'XTickLabel',[xticklab xticklab]);
        end
        if doCluster==1
            nClust=2;
            takeRows=size(tempSpectra',1);
            unitByUnitResponse=nan(length(psth.psths),takeRows*size(concatSpectra,2));
            refToUnits=nan(size(concatSpectra,1),1);
            for q=1:length(psth.psths)
                takeUnitSet=concatSpectra((q-1)*takeRows+1:(q-1)*takeRows+takeRows,:);
                takeUnitSet=takeUnitSet(~any(isnan(takeUnitSet),2));
                unitByUnitResponse(q,1:length(takeUnitSet(1:end)))=takeUnitSet(1:end);
                refToUnits((q-1)*takeRows+1:(q-1)*takeRows+takeRows)=q;
            end
            idx=kmeans(unitByUnitResponse(:,~any(isnan(unitByUnitResponse),1)),nClust);
            figure();
            new_concatSpectra=nan(size(concatSpectra));
            runningRow=0;
            allGrpAvs=cell(1,nClust);
            for q=1:nClust
                unitsInClust=find(idx==q);
                disp([repmat(q,length(unitsInClust),1) unitsInClust]);
                new_concatSpectra(runningRow+1:runningRow+length(unitsInClust)*size(tempSpectra',1),:)=concatSpectra(ismember(refToUnits,unitsInClust),:);
                currGrpAv=zeros(size([tempSpectra' tempSpectra2']));
                for r=1:length(unitsInClust)
                    curruic=unitsInClust(r);
                    currGrpAv=currGrpAv+concatSpectra(ismember(refToUnits,curruic),:);
                end
                currGrpAv=currGrpAv./length(unitsInClust);
                allGrpAvs{q}=currGrpAv;
                runningRow=runningRow+length(unitsInClust)*size(tempSpectra',1);
            end
            heatmap(new_concatSpectra);  
            figure();
            togAvs=[];
            for q=1:nClust
                togAvs=[togAvs; allGrpAvs{q}];
            end
            heatmap(togAvs);                
        end
    end
end
    
    

% for i=1:length(psth.psths)
%     temppre=allS{i};
% %     sumSpecgram=zeros(size(temppre{1}));
% %     for z=1:length(tryS)
% %         sumSpecgram=sumSpecgram+temppre{z};
% %     end
% %     sumSpecgram=sumSpecgram./length(tryS);
% %     currUnitMean=nanmean(sumSpecgram,1);
%     for z=1:length(tryS)
%         temp=temppre{z}-repmat(currUnitMean{i},size(temppre{z},1),1);
%         for q=1:size(temp,1)
%             temp(q,:)=temp(q,:)-min(temp(q,:));
%             temp(q,:)=temp(q,:)./max(temp(q,:));
%         end
%         out{z}=temp;
%     end
%     allS{i}=out;
% end

lowBand=[2.5 12];
highBand=[12 18];
allratios=cell(1,length(psth.psths));
for i=1:length(psth.psths)
    for z=1:length(tryS)
        temppre=allS{i};
        temp=temppre{z};
        currratio(z,:)=nanmean(temp(:,finterp>=highBand(1) & finterp<=highBand(2)),2)./nanmean(temp(:,finterp>=lowBand(1) & finterp(lowBand(2))),2);
    end
    allratios{i}=currratio;
end

figure();
% xticklab={'250','500','750','1000','1250','1500'};
for z=1:length(tryS)
    subplot(length(tryS),1,z);
    togetherRatios=[];
    for i=1:length(allratios)
        currratio=allratios{i};
        togetherRatios(i,:)=currratio(z,:);
    end
    nClust=4;
    idx=kmeans(togetherRatios,nClust);
    new_heatmap=[];
    for i=1:nClust
        new_heatmap=[new_heatmap; togetherRatios(idx==i,:)];
    end
    heatmap(new_heatmap); 
end
% h=gca;
% set(h,'XTick',[1:6]);
% set(h,'XTickLabel',xticklab);
% xlabel('Time (ms)');

end

function [f,t,S]=scriptForPlottingAutocorr_sub(psth,autocorr_out,usel,uses,i,lagsWindow,fromPSTH)

l=psth.unitLED{1};
s=psth.unitStimcond{1};
useTrials1=[];
for j=1:length(usel)
    useTrials1=[useTrials1 find(l>=usel(j)-0.001 & l<=usel(j)+0.001)];
end
useTrials1=sort(useTrials1);
useTrials2=[];
for j=1:length(uses)
    useTrials2=[useTrials2 find(s>=uses(j)-0.001 & s<=uses(j)+0.001)];
end
useTrials2=sort(useTrials2);
useTrials=useTrials1(ismember(useTrials1,useTrials2));
 
temp=psth.psths{i};

% figure(); 
temp=autocorr_out.autocorrs{i};
a=nanmean(temp(useTrials,:),1);
% plot(autocorr_out.lags.*20,a,'Color','k');
% ylim([0 max(a(autocorr_out.lags.*20>0))]);

params.Fs=50;
% params.tapers=[5 9];
params.tapers=[2 5];
params.trialave=1;
% params.fpass=[2.5 20];
params.fpass=[2.5 30];
temp(temp<10^-5)=0;
% movingwin=[0.3 0.02];
movingwin=[1 0.01];
% [S,f]=mtspectrumpb(nanmean(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2)),1),params);
% if any(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2))>0)
%     [S,f]=mtspectrumpb(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2)),params);
    if fromPSTH==1
        temp=psth.psths{i};
        [S,t,f]=mtspecgrampb([zeros(size(temp(useTrials,psth.t>=lagsWindow(1) & psth.t<=lagsWindow(2)),1),25) temp(useTrials,psth.t>=lagsWindow(1) & psth.t<=lagsWindow(2)) zeros(size(temp(useTrials,psth.t>=lagsWindow(1) & psth.t<=lagsWindow(2)),1),25)]',movingwin,params);
%         [S,t,f]=mtspecgrampb(temp(useTrials,psth.t>=lagsWindow(1) & psth.t<=lagsWindow(2))',movingwin,params);
    else
        [S,t,f]=mtspecgrampb(temp(useTrials,autocorr_out.lags.*20>lagsWindow(1) & autocorr_out.lags.*20<lagsWindow(2))',movingwin,params);
    end
% else
%     f=[params.fpass(1):0.5:params.fpass(2)];
%     S=zeros(size(f));
% end
% figure(); plot(f,S,'g');
% disp('hi')
end