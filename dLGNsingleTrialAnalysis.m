function dLGNsingleTrialAnalysis(datadir,workspace)

takeLED=true;
maxFreq=330; % in Hz
maxTime=12; % in seconds
maxFreq_override=50; % in Hz
maxTime_override=12; % in seconds from start
alphaRange=[10 14];
% alphaRange=[20 40];
% alphaRange=[5 8];
% sortingWindow=[3 3.6]; % time window used to sort trials into low and high
% sortingWindow=[5.5 6.5]; % time window used to sort trials into low and high
sortingWindow=[4 5]; % time window used to sort trials into low and high
visEvWindow=[4.5 6];

if isempty(workspace)
    if iscell(datadir)
        
        all_units_noTheta_lowF1=[];
        all_units_theta_lowF1=[];
        all_units_noTheta_LED_lowF1=[];
        all_units_theta_LED_lowF1=[];
        
        all_units_noTheta_highF1=[];
        all_units_theta_highF1=[];
        all_units_noTheta_LED_highF1=[];
        all_units_theta_LED_highF1=[];
        
        all_alpha_noTheta={};
        all_meanRate_noTheta={};
        all_visF1_noTheta={};
        all_alpha_theta={};
        all_meanRate_theta={};
        all_visF1_theta={};
        
        all_alpha_noTheta_LED={};
        all_meanRate_noTheta_LED={};
        all_visF1_noTheta_LED={};
        all_alpha_theta_LED={};
        all_meanRate_theta_LED={};
        all_visF1_theta_LED={};
               
        for i=1:length(datadir)
            d=datadir{i};
            disp(d);
           
            a=load([d '\noTheta_noLED.mat']);
            noTheta=a.noTheta;
            a=load([d '\theta_noLED.mat']);
            theta=a.theta;
            
            % collapseAcrossTrials
            [lowAlphaSpecs,highAlphaSpecs,alpha_noTheta,meanRate_noTheta,visF1_noTheta]=splitTrialsbyAlpha(noTheta,alphaRange,sortingWindow,sortingWindow,visEvWindow);
            all_alpha_noTheta(length(all_alpha_noTheta)+1:length(all_alpha_noTheta)+length(alpha_noTheta))=alpha_noTheta(1:end);
            all_meanRate_noTheta(length(all_meanRate_noTheta)+1:length(all_meanRate_noTheta)+length(meanRate_noTheta))=meanRate_noTheta(1:end);
            all_visF1_noTheta(length(all_visF1_noTheta)+1:length(all_visF1_noTheta)+length(visF1_noTheta))=visF1_noTheta(1:end);
            [lowAlphaSpecs_theta,highAlphaSpecs_theta,alpha_theta,meanRate_theta,visF1_theta]=splitTrialsbyAlpha(theta,alphaRange,sortingWindow,sortingWindow,visEvWindow);
            all_alpha_theta(length(all_alpha_theta)+1:length(all_alpha_theta)+length(alpha_theta))=alpha_theta(1:end);
            all_meanRate_theta(length(all_meanRate_theta)+1:length(all_meanRate_theta)+length(meanRate_theta))=meanRate_theta(1:end);
            all_visF1_theta(length(all_visF1_theta)+1:length(all_visF1_theta)+length(visF1_theta))=visF1_theta(1:end);
            
            if takeLED==true
                a=load([d '\noTheta_LED.mat']);
                noTheta_LED=a.noTheta;
                a=load([d '\theta_LED.mat']);
                theta_LED=a.theta;
                
                % collapseAcrossTrials
                [lowAlphaSpecs_LED,highAlphaSpecs_LED,alpha_noTheta,meanRate_noTheta,visF1_noTheta]=splitTrialsbyAlpha(noTheta,alphaRange,sortingWindow,sortingWindow,visEvWindow);
                all_alpha_noTheta_LED(length(all_alpha_noTheta_LED)+1:length(all_alpha_noTheta_LED)+length(alpha_noTheta))=alpha_noTheta(1:end);
                all_meanRate_noTheta_LED(length(all_meanRate_noTheta_LED)+1:length(all_meanRate_noTheta_LED)+length(meanRate_noTheta))=meanRate_noTheta(1:end);
                all_visF1_noTheta_LED(length(all_visF1_noTheta_LED)+1:length(all_visF1_noTheta_LED)+length(visF1_noTheta))=visF1_noTheta(1:end);
                [lowAlphaSpecs_theta_LED,highAlphaSpecs_theta_LED,alpha_theta,meanRate_theta,visF1_theta]=splitTrialsbyAlpha(theta,alphaRange,sortingWindow,sortingWindow,visEvWindow);
                all_alpha_theta_LED(length(all_alpha_theta_LED)+1:length(all_alpha_theta_LED)+length(alpha_theta))=alpha_theta(1:end);
                all_meanRate_theta_LED(length(all_meanRate_theta_LED)+1:length(all_meanRate_theta_LED)+length(meanRate_theta))=meanRate_theta(1:end);
                all_visF1_theta_LED(length(all_visF1_theta_LED)+1:length(all_visF1_theta_LED)+length(visF1_theta))=visF1_theta(1:end);
            end

            all_units_noTheta_lowF1=cat(3,all_units_noTheta_lowF1,lowAlphaSpecs);
            all_units_theta_lowF1=cat(3,all_units_theta_lowF1,lowAlphaSpecs_theta);
            all_units_noTheta_highF1=cat(3,all_units_noTheta_highF1,highAlphaSpecs);
            all_units_theta_highF1=cat(3,all_units_theta_highF1,highAlphaSpecs_theta);
            
            if takeLED==true
                all_units_noTheta_LED_lowF1=cat(3,all_units_noTheta_LED_lowF1,lowAlphaSpecs_LED);
                all_units_theta_LED_lowF1=cat(3,all_units_theta_LED_lowF1,lowAlphaSpecs_theta_LED);
                all_units_noTheta_LED_highF1=cat(3,all_units_noTheta_LED_highF1,highAlphaSpecs_LED);
                all_units_theta_LED_highF1=cat(3,all_units_theta_LED_highF1,highAlphaSpecs_theta_LED);
            end
        end
    else
        error('expected datadir to be a cell array');
    end
else
    load(workspace);
    maxFreq=maxFreq_override;
    maxTime=maxTime_override;
    all_units_noTheta_lowF1=backup_all_units_noTheta_lowF1;
    all_units_theta_lowF1=backup_all_units_theta_lowF1;
    all_units_noTheta_highF1=backup_all_units_noTheta_highF1;
    all_units_theta_highF1=backup_all_units_theta_highF1;
    all_units_noTheta_LED_lowF1=backup_all_units_noTheta_LED_lowF1;
    all_units_theta_LED_lowF1=backup_all_units_theta_LED_lowF1;
    all_units_noTheta_LED_highF1=backup_all_units_noTheta_LED_highF1;
    all_units_theta_LED_highF1=backup_all_units_theta_LED_highF1;
end

backup_all_units_noTheta_lowF1=all_units_noTheta_lowF1;
backup_all_units_theta_lowF1=all_units_theta_lowF1;
backup_all_units_noTheta_highF1=all_units_noTheta_highF1;
backup_all_units_theta_highF1=all_units_theta_highF1;

backup_all_units_noTheta_LED_lowF1=all_units_noTheta_LED_lowF1;
backup_all_units_theta_LED_lowF1=all_units_theta_LED_lowF1;
backup_all_units_noTheta_LED_highF1=all_units_noTheta_LED_highF1;
backup_all_units_theta_LED_highF1=all_units_theta_LED_highF1;

% Do whitening

% Norm

t=noTheta.allS.t;
f=noTheta.allS.f;

% Get alpha power vs F1 vs visual response, etc.
diffThresh=0; 
[useUnits_vis_noTheta,diffEvMinusSpont_noTheta]=isVisuallyResponsive(all_units_noTheta_lowF1,all_units_noTheta_highF1,[3 3.6],[4 6.5],t,f,diffThresh,[2.5 3.5]);
% useUnits_vis_noTheta=diffEvMinusSpont_noTheta<=1.25 & diffEvMinusSpont_noTheta>0;

versusRange=[20 40];
[alphaDiff,F1responseForGLM,meanRateForGLM,visEvForGLM,varF1ForGLM]=plotXAsFuncOfY(all_units_noTheta_lowF1,all_units_noTheta_highF1,alphaRange,versusRange,sortingWindow,noTheta,maxTime,[-0.8+0.01:0.02:0.8-0.01],[0 50],diffEvMinusSpont_noTheta);

[meanRateDiff,F1response]=plotXAsFuncOfY(all_units_noTheta_lowF1,all_units_noTheta_highF1,[0 50],[],sortingWindow,noTheta,maxTime,[-0.8+0.01:0.02:0.8-0.01],[],[]);

% mdl=fitglm([alphaDiff meanRateForGLM visEvForGLM],F1responseForGLM,'interactions')

mdl=fitglm([alphaDiff meanRateForGLM visEvForGLM],varF1ForGLM,'interactions')

figure(); scatter(alphaDiff,F1responseForGLM);
ranksum(F1responseForGLM(visEvForGLM>10),F1responseForGLM(visEvForGLM<10))

%takeUnits=1:765;
takeUnits=find(useUnits_vis_noTheta);
temp1=all_units_noTheta_lowF1;
for i=1:6
    temp1(:,i,:)=temp1(:,i,:)-repmat(nanmean(temp1(t>=0 & t<=3,i,:),1),size(temp1,1),1,1);
    temp1(:,i,:)=temp1(:,i,:)+repmat(nanmean(nanmean(temp1(t>=0 & t<=3,7:end,:),1),2),size(temp1,1),1,1)-0.2;
end
temp2=all_units_noTheta_highF1;
for i=1:6
    temp2(:,i,:)=temp2(:,i,:)-repmat(nanmean(temp2(t>=0 & t<=3,i,:),1),size(temp2,1),1,1);
    temp2(:,i,:)=temp2(:,i,:)+repmat(nanmean(nanmean(temp2(t>=0 & t<=3,7:end,:),1),2),size(temp2,1),1,1)-0.2;
end
plotSpecgramsComparingHighAndLow(temp1,temp2,takeUnits,t,f,maxTime,maxFreq,'no theta');
powerInBand_low=takeFrequencyBand(temp1(t<=maxTime,:,takeUnits),[2.5 3.5],noTheta);
powerInBand_high=takeFrequencyBand(temp2(t<=maxTime,:,takeUnits),[2.5 3.5],noTheta);
plotWStderr(t(t<=maxTime),powerInBand_low,powerInBand_high,'k','r',0);
[alphaDiff,F1responseForGLM,meanRateForGLM,visEvForGLM,varF1ForGLM]=plotXAsFuncOfY(all_units_noTheta_lowF1(:,:,takeUnits),all_units_noTheta_highF1(:,:,takeUnits),alphaRange,versusRange,sortingWindow,noTheta,maxTime,[-0.8+0.01:0.02:0.8-0.01],[0 50],diffEvMinusSpont_noTheta(takeUnits));
mdl=fitglm([alphaDiff meanRateForGLM visEvForGLM],F1responseForGLM,'interactions')

takeUnits=187:540;
plotSpecgramsComparingHighAndLow(all_units_noTheta_lowF1-0.229,all_units_noTheta_highF1,takeUnits,t,f,maxTime,maxFreq,'no theta');
powerInBand_low=takeFrequencyBand(all_units_noTheta_lowF1(t<=maxTime,:,takeUnits)-0.229,[2.5 3.5],noTheta);
powerInBand_high=takeFrequencyBand(all_units_noTheta_highF1(t<=maxTime,:,takeUnits),[2.5 3.5],noTheta);
plotWStderr(t(t<=maxTime),powerInBand_low,powerInBand_high,'k','r',0);
signrank(nanmean(powerInBand_low(:,t(t<=maxTime)>4 & t(t<=maxTime)<=6),2)-nanmean(powerInBand_low(:,t(t<=maxTime)>2.5 & t(t<=maxTime)<=3.7),2),nanmean(powerInBand_high(:,t(t<=maxTime)>4 & t(t<=maxTime)<=6),2)-nanmean(powerInBand_high(:,t(t<=maxTime)>2.5 & t(t<=maxTime)<=3.7),2))

end

function [lowAlphaSpecs,highAlphaSpecs,alpha_noTheta,meanRate_noTheta,visF1_noTheta]=splitTrialsbyAlpha(noTheta,alphaRange,alphaSortingWindow,meanSortingWindow,visSortingWindow)

takePercentiles=true;

for i=1:length(noTheta.allS.S)
    temp=noTheta.allS.S{i};
    [alpha_noTheta{i},meanRate_noTheta{i},visF1_noTheta{i}]=getPowerInBands(temp,alphaRange,noTheta,alphaSortingWindow,meanSortingWindow,visSortingWindow); % dim 3 is trials
    if takePercentiles==false
        media=nanmedian(alpha_noTheta{i}); % dim 1 is now trials
        takeHigh=alpha_noTheta{i}>=media;
        takeLow=alpha_noTheta{i}<media;
    else
%         takeHigh=alpha_noTheta{i}>prctile(alpha_noTheta{i},65);
%         takeLow=alpha_noTheta{i}<prctile(alpha_noTheta{i},35);
        takeHigh=alpha_noTheta{i}>prctile(alpha_noTheta{i},80);
        takeLow=alpha_noTheta{i}<prctile(alpha_noTheta{i},20);
    end
    lowAlphaSpecs(:,:,i)=reshape(nanmean(temp(:,:,takeLow),3),size(temp,1),size(temp,2));
    highAlphaSpecs(:,:,i)=reshape(nanmean(temp(:,:,takeHigh),3),size(temp,1),size(temp,2));
end

end

function [alpha_noTheta,meanRate_noTheta,visF1_noTheta]=getPowerInBands(all_units_noTheta,alphaRange,noTheta,alphaTime,meanTime,visTime)

versusRange=[20 40];
% versusRange=[0 50];
alpha_noTheta=takeFrequencyBand(all_units_noTheta,alphaRange,noTheta);
alpha_noTheta=nanmean(alpha_noTheta(:,noTheta.allS.t>=alphaTime(1) & noTheta.allS.t<=alphaTime(2)),2); % dim 1 is trials
versus_noTheta=takeFrequencyBand(all_units_noTheta,versusRange,noTheta);
versus_noTheta=nanmean(versus_noTheta(:,noTheta.allS.t>=alphaTime(1) & noTheta.allS.t<=alphaTime(2)),2);
alpha_noTheta=alpha_noTheta./versus_noTheta;

meanRate_noTheta=takeFrequencyBand(all_units_noTheta,[0 50],noTheta);
meanRate_noTheta=nanmean(meanRate_noTheta(:,noTheta.allS.t>=meanTime(1) & noTheta.allS.t<=meanTime(2)),2);

visF1_noTheta=takeFrequencyBand(all_units_noTheta,[2.5 3.5],noTheta);
visF1_noTheta=nanmean(visF1_noTheta(:,noTheta.allS.t>=visTime(1) & noTheta.allS.t<=visTime(2)),2);

end

function [concatX,concatY,concatThird,concatParam,concatVar]=plotXAsFuncOfY(all_units_low,all_units_high,useFreqRange,versusRange,useTimeRange,noTheta,maxTime,tryThresh,returnThirdVarRange,param)

useUnits_noTheta_alphaMatch=cell(1,length(tryThresh)-1);
for i=1:length(tryThresh)-1
    diffThresh=tryThresh(i);
    [~,diffHighMinusLow_noTheta]=isAlphaHigherForHigh(all_units_low,all_units_high,useFreqRange,useTimeRange,noTheta.allS.t,noTheta.allS.f,diffThresh,versusRange);
    useUnits_noTheta_alphaMatch{i}=diffHighMinusLow_noTheta>=tryThresh(i) & diffHighMinusLow_noTheta<=tryThresh(i+1);
end
   
t=noTheta.allS.t;
powLow=cell(1,length(tryThresh)-1);
powHigh=cell(1,length(tryThresh)-1);
powThird=cell(1,length(tryThresh)-1);
highMinusLow=cell(1,length(tryThresh)-1);
powParam=cell(1,length(tryThresh)-1);
for i=1:length(tryThresh)-1
    useUnits_noTheta=useUnits_noTheta_alphaMatch{i};
    powerInBand_low=takeFrequencyBand(all_units_low(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta);
    powLow{i}=nanmean(powerInBand_low(:,t>=4 & t<=6.5),2);
    powerInBand_high=takeFrequencyBand(all_units_high(t<=maxTime,:,useUnits_noTheta),[2.5 3.5],noTheta);
    powHigh{i}=nanmean(powerInBand_high(:,t>=4 & t<=6.5),2);
    highMinusLow{i}=powHigh{i}-powLow{i};
    if ~isempty(returnThirdVarRange)
        powerInBand_third_low=takeFrequencyBand(all_units_low(t<=maxTime,:,useUnits_noTheta),returnThirdVarRange,noTheta);
        powerInBand_third_high=takeFrequencyBand(all_units_high(t<=maxTime,:,useUnits_noTheta),returnThirdVarRange,noTheta);
        powThird{i}=nanmean(powerInBand_third_high(:,t>=4 & t<=6.5),2)-nanmean(powerInBand_third_low(:,t>=4 & t<=6.5),2);
    end
    if ~isempty(param)
        powParam{i}=param(useUnits_noTheta)';
    end
end

figure();
concatX=[];
concatY=[];
concatThird=[];
concatParam=[];
isInParamRange=[];
concatVar=[];
for i=1:length(tryThresh)-1
    if ~isempty(param)
        paramThresh=0;
        paramThresh2=1.25;
        paramTemp=powParam{i};
        takeParaming=paramTemp>paramThresh & paramTemp<=paramThresh2;
        x=ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2);
        y=highMinusLow{i};
        scatter(x(~takeParaming),y(~takeParaming),[],'k'); hold all;
        scatter(x(takeParaming),y(takeParaming),[],'r');
        scatter(x(paramTemp>paramThresh2),y(paramTemp>paramThresh2),[],'g');
    else
        scatter(ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2),highMinusLow{i});
    end
    concatX=[concatX; ones(size(highMinusLow{i})).*nanmean([tryThresh(i) tryThresh(i+1)],2)];
    concatY=[concatY; highMinusLow{i}];
    concatThird=[concatThird; powThird{i}];
    if ~isempty(param)
        concatParam=[concatParam; powParam{i}];
        isInParamRange=[isInParamRange; takeParaming];
    end
    concatVar=[concatVar; nanvar(highMinusLow{i},1).*ones(size(highMinusLow{i}))];
    hold on;
end

[r,p]=corrcoef(concatX,concatY)

[n,x]=hist(concatY(isInParamRange==1),100);
figure(); plot(x,n./nansum(n),'Color','r');
hold on;
[n,x]=hist(concatY(isInParamRange==0),300);
plot(x,n./nansum(n),'Color','k');

end

function [useUnits,diffEvMinusSpont]=isVisuallyResponsive(all_units_low,all_units_high,spontWindow,evWindow,t,f,diffThresh,visRange)

if length(size(all_units_low))>2
    spont=0.5*reshape(nanmean(nanmean(all_units_low(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_low,3)) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_high,3));
    ev=0.5*reshape(nanmean(nanmean(all_units_low(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_low,3)) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2),:),2),1),1,size(all_units_high,3));
else
    spont=0.5*reshape(nanmean(nanmean(all_units_low(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=spontWindow(1) & t<=spontWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1);
    ev=0.5*reshape(nanmean(nanmean(all_units_low(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1) + 0.5*reshape(nanmean(nanmean(all_units_high(t>=evWindow(1) & t<=evWindow(2),f>=visRange(1) & f<=visRange(2)),2),1),1,1);
end
if isempty(diffThresh)
    useUnits=spont<ev;
    diffEvMinusSpont=ev-spont;
else
    diffEvMinusSpont=ev-spont;
    useUnits=diffEvMinusSpont>diffThresh;
end


end

function plotSpecgramsComparingHighAndLow(all_units_noTheta_lowF1_Ntsr1,all_units_noTheta_highF1_Ntsr1,useUnits_noTheta,t,f,maxTime,maxFreq,titAdd)

f=f(2:end);

test=reshape(nanmean(all_units_noTheta_lowF1_Ntsr1(:,2:end,useUnits_noTheta),3),size(all_units_noTheta_lowF1_Ntsr1,1),size(all_units_noTheta_lowF1_Ntsr1,2)-1);
test_high=reshape(nanmean(all_units_noTheta_highF1_Ntsr1(:,2:end,useUnits_noTheta),3),size(all_units_noTheta_highF1_Ntsr1,1),size(all_units_noTheta_highF1_Ntsr1,2)-1);
figure(); imagesc(t(t<maxTime),f(f<=maxFreq),test_high(t<maxTime,f<=maxFreq)'-test(t<maxTime,f<=maxFreq)'); title(['high minus low ' titAdd]);
figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test_high(t<maxTime,f<=maxFreq) test(t<maxTime,f<=maxFreq)]'); title(['high then low ' titAdd]);
figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test_high(t<maxTime,f<=maxFreq)]'); title(['high ' titAdd]);
figure(); imagesc(t(t<maxTime),f(f<=maxFreq),[test(t<maxTime,f<=maxFreq)]'); title(['low ' titAdd]);

end

function [useUnits,diffHighMinusLow]=isAlphaHigherForHigh(all_units_low,all_units_high,alphaRange,timeRange,t,f,diffThresh,versusRange)

if length(size(all_units_low))>2
    alphaForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2),:),2),1),1,size(all_units_low,3));
    alphaForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2),:),2),1),1,size(all_units_high,3));
else
    alphaForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2)),2),1),1,1);
    alphaForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=alphaRange(1) & f<=alphaRange(2)),2),1),1,1);
end
if ~isempty(versusRange)
    if length(size(all_units_low))>2
        versusForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2),:),2),1),1,size(all_units_low,3));
        versusForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2),:),2),1),1,size(all_units_high,3));
    else
        versusForLow=reshape(nanmean(nanmean(all_units_low(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2)),2),1),1,1);
        versusForHigh=reshape(nanmean(nanmean(all_units_high(t>=timeRange(1) & t<=timeRange(2),f>=versusRange(1) & f<=versusRange(2)),2),1),1,1);
    end
    alphaForLow=alphaForLow./versusForLow;
    alphaForHigh=alphaForHigh./versusForHigh;
end
if isempty(diffThresh)
    useUnits=alphaForLow<alphaForHigh;
    diffHighMinusLow=alphaForHigh-alphaForLow;
else
    diffHighMinusLow=alphaForHigh-alphaForLow;
    useUnits=diffHighMinusLow>diffThresh;
end
    
end

function [y1,y2]=plotWStderr(x,y1,y2,c1,c2,doFill)

baseSub=false;
baseWindow=[1.485 3.535];

if baseSub==true
    y1=y1-repmat(nanmean(y1(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y1,2));
    if ~isempty(y2)
        y2=y2-repmat(nanmean(y2(:,x>=baseWindow(1) & x<=baseWindow(2)),2),1,size(y2,2));
    end
end

figure();
plot(x,nanmean(y1,1),'Color',c1);
hold on;
if doFill==1
    fill([x fliplr(x)],[nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)) fliplr(nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)))],[0.5 0.5 0.5]);
end
plot(x,nanmean(y1,1),'Color',c1);
plot(x,nanmean(y1,1)+nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);
plot(x,nanmean(y1,1)-nanstd(y1,[],1)./sqrt(size(y1,1)),'Color',c1);

if ~isempty(y2)
    plot(x,nanmean(y2,1),'Color',c2);
    if doFill==1
        fill([x fliplr(x)],[nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)) fliplr(nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)))],[0.1 0.7 0.5]);
    end
    plot(x,nanmean(y2,1),'Color',c2);
    plot(x,nanmean(y2,1)+nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    plot(x,nanmean(y2,1)-nanstd(y2,[],1)./sqrt(size(y2,1)),'Color',c2);
    
    plot(x,nanmean(y1,1),'Color',c1);
end
    
end

function powerInBand=takeTimeSlice(specgram,timeWindow,noTheta)

if length(size(specgram))>2
    powerInBand=nan(size(specgram,3),size(specgram,2)); 
    for i=1:size(specgram,3)
        powerInBand(i,:)=reshape(nanmean(specgram(noTheta.allS.t>=timeWindow(1) & noTheta.allS.t<=timeWindow(2),:,i),1),size(specgram,2),1);
    end
else
    powerInBand=nan(1,size(specgram,2)); 
    powerInBand(1,:)=reshape(nanmean(specgram(noTheta.allS.t>=timeWindow(1) & noTheta.allS.t<=timeWindow(2),:),1),size(specgram,2),1);
end

end

function powerInBand=takeFrequencyBand(specgram,freqBand,noTheta)

if length(size(specgram))>2
    powerInBand=nan(size(specgram,3),size(specgram,1)); 
    for i=1:size(specgram,3)
        powerInBand(i,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2),i),2),size(specgram,1),1);
    end
else
    powerInBand=nan(1,size(specgram,1)); 
    powerInBand(1,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2)),2),size(specgram,1),1);
end

end

function combined=addSigMask(data,sigMask,noTheta,p)

tempfill=nanmean(nanmean(data,1),2);
combined=data.*(sigMask(noTheta.allS.t<12,noTheta.allS.f<=50)<p);
combined(combined==0)=tempfill;

end

function sigMask=calculateSignificanceMask(dataset1,dataset2)

% dim 3 is units

sigMask=nan(size(dataset1,1),size(dataset2,2));
for i=1:size(dataset1,1)
    for j=1:size(dataset1,2)
        temp1=dataset1(i,j,:);
        temp2=dataset2(i,j,:);
        if all(isnan(temp1)) || all(isnan(temp2))
            sigMask(i,j)=1;
        else
            sigMask(i,j)=signrank(reshape(temp1,size(temp1,3),1),reshape(temp2,size(temp2,3),1));        
        end
    end
end

end

function output=normAll(input)

for i=1:size(input,3)
    temp=reshape(input(:,:,i),size(input,1),size(input,2));
    normTo=nansum(temp(:,5));
    temp(:,1)=(temp(:,1)./nansum(temp(:,1))).*normTo;
    temp(:,2)=(temp(:,2)./nansum(temp(:,2))).*normTo;
    temp(:,3)=(temp(:,3)./nansum(temp(:,3))).*normTo;
    temp(:,4)=(temp(:,4)./nansum(temp(:,4))).*normTo;
    temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    output(:,:,i)=temp;
end

end

function output=whitenAll(input,isCx,isHigh)

if isCx==false
    for i=1:size(input,3)
        temp=reshape(input(:,:,i),size(input,1),size(input,2));
        temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        temp(:,3)=1.0475*temp(:,3);
        output(:,:,i)=temp;
    end
else
    for i=1:size(input,3)
%         temp=reshape(input(:,:,i),size(input,1),size(input,2));
%         temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*0.37);
%         temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.39);
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,1)=1.05*temp(:,1);
%         temp(:,3)=1.06*temp(:,3);
%         temp(:,4)=0.75*temp(:,4);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
%         output(:,:,i)=temp;
         if isHigh==0
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.1);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.1);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.48*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
            
         else
%             temp=reshape(input(:,:,i),size(input,1),size(input,2));
%             temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1.5);
%             temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.9);
%             temp(:,2)=1.07*temp(:,2);
%             temp(:,1)=1.03*temp(:,1);
%             temp(:,3)=1.06*temp(:,3);
%             temp(:,4)=0.75*temp(:,4);
%             output(:,:,i)=temp;
            
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.3);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.3);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.43*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
         end
    end
end

end

function temp=whitenAndNormSpecgram(temp,isCx,isHigh)

if isCx==false
    temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
    temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
    temp(:,2)=1.07*temp(:,2);
    temp(:,4)=1.07*temp(:,4);
    temp(:,3)=1.0475*temp(:,3);
    % temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    % temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
else
    if isHigh==0
%         temp(:,1:3)=temp(:,1:3)./(3.3*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.0*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));


%         % BEST FOR NO THETA
%         temp(:,1:3)=temp(:,1:3)./(3.53*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.13*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));


        temp(:,1:3)=temp(:,1:3)./(4*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(3*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        %     temp(:,3)=1.0475*temp(:,3);
        temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
        temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    else
%         temp(:,1:3)=temp(:,1:3)./(3.3*0.17.*(0.7225/0.1121)*1.2);
%         temp(:,4)=temp(:,4)./(2.0*0.73.*(0.7225/0.5035)*0.9);
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.03*temp(:,4);
%         temp(:,3)=1.0*temp(:,3);
%         temp(:,4)=0.75*temp(:,4);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));

%         % BEST FOR NO THETA
%         temp(:,1:3)=temp(:,1:3)./(3.53*0.17.*(0.7225/0.1121));
%         temp(:,4)=temp(:,4)./(2.13*0.73.*(0.7225/0.5035));
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,4)=1.07*temp(:,4);
%         %     temp(:,3)=1.0475*temp(:,3);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));

        temp(:,1:3)=temp(:,1:3)./(4*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(3*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        %     temp(:,3)=1.0475*temp(:,3);
        temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
        temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
    end
end

K=ones([2 2]); 
temp=conv2(temp,K,'same');
temp=interp2(temp,4);

temp=temp(1:end-13,1:end-13);

end

function out=addAll(data1,data2)

tmp=cat(3,data1,data2);
out=nansum(tmp,3);

end

function align_specgram=alignToStim(specgram,specgram_t,placeStimOnsetAt,currStimOnset)

if placeStimOnsetAt==currStimOnset
    align_specgram=specgram;
    return
end

timeDiff=abs(currStimOnset-placeStimOnsetAt);
inds_diff=floor(timeDiff/(specgram_t(2)-specgram_t(1)));

if currStimOnset>placeStimOnsetAt
    align_specgram=[specgram(inds_diff:end,:); nan(inds_diff-1,size(specgram,2))];
elseif currStimOnset<placeStimOnsetAt
    align_specgram=[nan(inds_diff,size(specgram,2)); specgram(1:end-inds_diff,:)];
end

end

function a=loadVar(varargin)

if length(varargin)==3
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=[];
elseif length(varargin)==4
    varDir=varargin{1};
    varName=varargin{2};
    loadFieldByField=varargin{3};
    inFileName=varargin{4};
else
    error('expected 3 or 4 arguments to loadVar');
end

if loadFieldByField==true
    if isempty(inFileName)
        a.(varName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(varName))
            a=a.(varName);
        end
    else
        a.(inFileName)=loadStructFieldByField([varDir '\' varName]);
        if iscell(a.(inFileName))
            a=a.(inFileName);
        end
    end
else
    a=load([varDir '\' varName '.mat']);
end

end
       