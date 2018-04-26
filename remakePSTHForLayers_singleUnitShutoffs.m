function [psthsAcrossLayers,allUnitsTimeCourse]=remakePSTHForLayers_singleUnitShutoffs(bestSpikes,newassignsinfo,exptData,params,saveDir)

exptDir='W:\Analysis Computer\141003 Expts\'; %'W:\Expt Backups\Experiments\';

uniTrodes=unique(newassignsinfo.trode);

% Load expt and get real channels for these trodes
betterAssignsInfo=newassignsinfo;
cload=load([exptDir exptData.name '_expt.mat']);
expt=cload.expt; 
modLessThanOne=logical(mod(betterAssignsInfo.calibrated_evCh,4)<1);
modGreaterThanOne=logical(mod(betterAssignsInfo.calibrated_evCh,4)>=1);
for i=1:length(uniTrodes)
    cUniTrode=uniTrodes(i);
    thisChs=expt.sort.trode(cUniTrode).channels;
    if thisChs(1)==23
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=1;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne),4)+4;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne),4)+0;
%         betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode),4)+0;
    elseif thisChs(1)==19
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=2;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne),4)+8;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne),4)+4;
%         betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode),4)+4;
    elseif thisChs(1)==17
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=3;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne),4)+12;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne),4)+8;
%         betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode),4)+8;
    elseif thisChs(1)==18
        betterAssignsInfo.trode(betterAssignsInfo.trode==cUniTrode)=4;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modLessThanOne),4)+16;
        betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode & modGreaterThanOne),4)+12;
%         betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode)=mod(betterAssignsInfo.calibrated_evCh(betterAssignsInfo.trode==cUniTrode),4)+12;
    else
        disp('problem with trode assigns');
    end
end

concatassigns=[];
newass_assigns=[];
uniTrodes=unique(betterAssignsInfo.trode);
maxSoFar=0;
for i=1:length(uniTrodes)
    currTrode=uniTrodes(i);
    if i==1
        concatassigns=betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode);
        newass_assigns=betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode);
    else
        maxSoFar=max(concatassigns);
        concatassigns=[concatassigns betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode)+maxSoFar];
        newass_assigns=[newass_assigns betterAssignsInfo.original_assigns(betterAssignsInfo.trode==currTrode)];
    end
end
disp('Should be the same numbers, unless bestUnitsInfo already fixed');
disp([unique(bestSpikes.assigns)' concatassigns']);
% disp(concatassigns);

% newcalev=[];
% for i=1:length(concatassigns)
%     getnewass_assigns=newass_assigns(i);
% %     getInd=find(newassignsinfo.original_assigns==getnewass_assigns);
% %     calee=newassignsinfo.calibrated_evCh(getInd);
%     calee=newassignsinfo.calibrated_evCh(newassignsinfo.original_assigns==getnewass_assigns);
%     newcalev=[newcalev calee];
% end
% betterAssignsInfo.calibrated_evCh=newcalev;

psthsAcrossLayers.y1=zeros(15,length(exptData.xpoints));
psthsAcrossLayers.y2=zeros(15,length(exptData.xpoints));
psthsAcrossLayers.x=exptData.xpoints;
for i=1:15
    if i==15
        useAss=concatassigns(betterAssignsInfo.calibrated_evCh>=i & betterAssignsInfo.calibrated_evCh<=i+1);
    else
        useAss=concatassigns(betterAssignsInfo.calibrated_evCh>=i & betterAssignsInfo.calibrated_evCh<i+1);
    end
    [x,y1,y2]=scriptForComparingLayers(filtspikes(bestSpikes,0,'assigns',useAss),exptData.AIs,[],[],exptData,params);
    psthsAcrossLayers.y1(i,:)=y1;
    psthsAcrossLayers.y2(i,:)=y2;
end

layerData=psthsAcrossLayers;
save([saveDir '\fromunits_layerDataBETTER.mat'],'layerData');

% Now calculate shut-off for single units
itAssigns=concatassigns;
allUnitsTimeCourse.x=NaN(1,length(x));
allUnitsTimeCourse.y1=NaN(length(itAssigns),length(x));
allUnitsTimeCourse.y2=NaN(length(itAssigns),length(x));
allUnitsTimeCourse.taus=NaN(length(itAssigns),1);
allUnitsTimeCourse.assigns=itAssigns;
allUnitsTimeCourse.calibrated_evCh=betterAssignsInfo.calibrated_evCh;
for i=1:length(itAssigns)
    useAss=itAssigns(i);
    unitTimeCourse=fitSingleUnitShutoff(bestSpikes,useAss,exptData,length(x),params);
    allUnitsTimeCourse.x=unitTimeCourse.x;
    allUnitsTimeCourse.y1(i,:)=unitTimeCourse.y1;
    allUnitsTimeCourse.y2(i,:)=unitTimeCourse.y2;
    allUnitsTimeCourse.taus(i)=unitTimeCourse.tau;
    allUnitsTimeCourse.assigns(i)=useAss;
end
save([saveDir '\units_shutoffBETTER.mat'],'allUnitsTimeCourse');

end

function unitTimeCourse=fitSingleUnitShutoff(spikes,useAss,exptData,defaultXLength,params)
ledValue=exptData.useLEDcond{2}; % red
noLedValues=exptData.useLEDcond{1}; % black
waitS=0.003; % in s
finishS=waitS+0.06; % in s
Awindow=[-0.05 0]; % in s
zeroWindow=[0.05 0.1]; % in s

ledOnset=exptData.stimulusOn(1)+exptData.stimulusDuration/1000;
try
    [x,y1,y2]=scriptForComparingLayers(filtspikes(spikes,0,'assigns',useAss),exptData.AIs,[],[],exptData,params);
    unitTimeCourse.x=x;
    unitTimeCourse.y1=y1;
    unitTimeCourse.y2=y2;
catch
    unitTimeCourse.x=NaN(1,defaultXLength);
    unitTimeCourse.y1=NaN(1,defaultXLength);
    unitTimeCourse.y2=NaN(1,defaultXLength);
    unitTimeCourse.tau=nan;
    return
end
x=x';
y1=y1';
y2=y2';
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

function [xpoints1,ypoints1,ypoints2]=scriptForComparingLayers(spikes,useTheseFileInds,useBlackTrials,useRedTrials,exptData,params)
ledValue=exptData.useLEDcond{2}; % red
noLedValues=exptData.useLEDcond{1}; % black
binsize=1; % ms

stimconds=exptData.useStimcond;

nunits=1;
spikes=filtspikes(spikes,0,'stimcond',stimconds); 
if ~isempty(useTheseFileInds)
    temp=[];
    temp1=[];
    for i=1:length(ledValue)
        spikes=makeTempField(spikes,'led',ledValue(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_withLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
    temp=[];
    temp1=[];
    for i=1:length(noLedValues)
        spikes=makeTempField(spikes,'led',noLedValues(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_noLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
else
    temp=[];
    temp1=[];
    for i=1:length(ledValue)
        spikes=makeTempField(spikes,'led',ledValue(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_withLED=filtspikes(spikes,0,'temp',1);
    temp=[];
    temp1=[];
    for i=1:length(noLedValues)
        spikes=makeTempField(spikes,'led',noLedValues(i));
        temp(i,:)=spikes.temp;
        temp1(i,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    allUnits_noLED=filtspikes(spikes,0,'temp',1);
end

allSpiketimes_withLED=[];
allSpiketimes_noLED=[];

if ~isempty(useBlackTrials)
    allUnits_noLED=filtspikes(allUnits_noLED,0,'trials',sort(useBlackTrials));
end
if ~isempty(useRedTrials)
    allUnits_withLED=filtspikes(allUnits_withLED,0,'trials',sort(useRedTrials));
end
    
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_forLayers(allUnits_noLED,binsize,0,params.trialDuration);

[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psth_forLayers(allUnits_withLED,binsize,0,params.trialDuration);

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;

end

function [varargout] = psth_forLayers(spikes,binsize,bsmooth,duration)

if nargin < 2
    binsize = 50; 
end

if nargin < 3
    bsmooth = 0;
end

% Set duration and number of trials
% duration = 5;     % Maybe add checking for equal sweep durations?
if isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
%     if length(spikes.sweeps.trials)>=4*length(a)-5 && length(spikes.sweeps.trials)<=4*length(a)+5
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% numtrials=length(unique(spikes.trials));
% numtrials = length(spikes.sweeps.trials);
% numtrials=119;

% Set spiketimes
spiketimes = spikes.spiketimes;

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

nsForStdev=zeros(length(unique(spikes.trials)),size(n,2));
allTrials=unique(spikes.trials);
for i=1:length(allTrials)
    cspikes=filtspikes(spikes,0,'trials',allTrials(i));
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;

if all(isnan(n))
    n = 0;
end

% Compute center of bins
centers = edges + diff(edges(1:2))/2;

% Last point of n contains values falling on edge(end) -- usually zero
if bsmooth
% hPsth = line('XData',centers(1:end-1),'YData',smooth(n(1:end-1),3),...
%     'Parent',hAxes,'LineWidth',1.5);
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
%     hPsth = line('XData',centers(1:end-1),'YData',n(1:end-1),...
%     'Parent',hAxes,'LineWidth',1.5);
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

% Set default axes properties
% if sum(n) == 0
%     maxN = 0.1;
% else
%     maxN = max(n);
% end
% axis([0 duration 0 maxN])
% set(hAxes,'TickDir','out','FontSize',8)
% xlabel(hAxes,'seconds')
% ylabel(hAxes,'spikes/s')

% Outputs
% varargout{1} = hPsth;
% varargout{2} = hAxes;
varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),1);

% figure();
% plot(centers(1:end-1),n(1:end-1),'Color','black');
% hold on; 
% plot(centers(1:end-1),n(1:end-1)-1*std(nsForStdev(1:end-1),1),'Color','green');
% plot(centers(1:end-1),n(1:end-1)+1*std(nsForStdev(1:end-1),1),'Color','green');
end