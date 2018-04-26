function estimate_taus=deconvolutionDynamics(name,V1spikes,thalspikes,psths,saveDir)

test=0;
showFigs=0;
doAlign=0;
doFlip=0;
baseSubtract=0;
delayShift=1;
onlyStim=0;

if isempty(saveDir)
    saveDir='W:\Analysis Computer\Frequency Response With More Data 150120\Deconvolution Analysis\dLGN Only Anesth Test\';
end
    
params.useLEDcond=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
params.duration=4;
params.binsize=1;
params.spontWindow=[0 1];
params.stimWindow=[0 1];

saveDirName=[saveDir name];
if exist(saveDirName,'dir')==0
    st=mkdir(saveDirName);
    if ~st
        disp('Could not create save directory');
    end
end
save([saveDirName '\params.mat'],'params');

if isempty(psths)
    for i=1:length(params.useLEDcond)
        disp('calculating PSTHs');
        ledValue=params.useLEDcond(i);
        temp=[];
        temp1=[];
        for j=1:length(ledValue)
            V1spikes=makeTempField(V1spikes,'led',ledValue(j));
            temp(j,:)=V1spikes.temp;
            temp1(j,:)=V1spikes.sweeps.temp;
            V1spikes.temp=[];
            V1spikes.sweeps.temp=[];
        end
        V1spikes.temp=sum(temp,1)>=1;
        V1spikes.sweeps.temp=sum(temp1,1)>=1;
        currV1Spikes=filtspikes(V1spikes,0,'temp',1);
        temp=[];
        temp1=[];
        for j=1:length(ledValue)
            thalspikes=makeTempField(thalspikes,'led',ledValue(j));
            temp(j,:)=thalspikes.temp;
            temp1(j,:)=thalspikes.sweeps.temp;
            thalspikes.temp=[];
            thalspikes.sweeps.temp=[];
        end
        thalspikes.temp=sum(temp,1)>=1;
        thalspikes.sweeps.temp=sum(temp1,1)>=1;
        currThalSpikes=filtspikes(thalspikes,0,'temp',1);
        [~,~,~,xpoints1,ypoints1]=psth_wStdev(currV1Spikes,params.binsize,0,params.duration);
        [~,~,~,xpoints1thal,ypoints1thal]=psth_wStdev(currThalSpikes,params.binsize,0,params.duration);
        psthsV1.V1_x{i}=xpoints1;
        psthsV1.V1_y{i}=ypoints1;
        psthsThal.thal_x{i}=xpoints1thal;
        psthsThal.thal_y{i}=ypoints1thal;
        psths.V1_x{i}=xpoints1;
        psths.V1_y{i}=ypoints1;
        psths.thal_x{i}=xpoints1thal;
        psths.thal_y{i}=ypoints1thal;
    end
else
    psthsV1.V1_x=psths.V1_x;
    psthsV1.V1_y=psths.V1_y;
    psthsThal.thal_x=psths.thal_x;
    psthsThal.thal_y=psths.thal_y;
end
save([saveDirName '\psths.mat'],'psths');
save([saveDirName '\psthsV1.mat'],'psthsV1');
save([saveDirName '\psthsThal.mat'],'psthsThal');

if doFlip==1
    cxx=psths.V1_x{1};
    for i=1:length(psths.V1_y)
    	cxy=psths.V1_y{i};
        thaly=psths.thal_y{i};
        [cc,cc_lags]=xcorr(cxy-mean(cxy),thaly-mean(thaly),[],'none');
        [~,maxcc]=max(cc);
        ind=cc_lags(maxcc);
%         if showFigs==1
%             figure();
%             plot(downSampAv(cxx,10),downSampAv(cxy,10),'Color','k');
%             hold on;
%             plot(downSampAv(cxx,10),downSampAv(thaly,10),'Color','b');
%         end
        if ind<0
            % Flip yes
            cxy=-cxy;
            psths.V1_y{i}=cxy;
        end
%         if showFigs==1
%             figure();
%             plot(downSampAv(cxx,10),downSampAv(cxy,10),'Color','r');
%             hold on;
%             plot(downSampAv(cxx,10),downSampAv(thaly,10),'Color','b');
%         end
    end
end
if baseSubtract==1
    cxx=psths.V1_x{1};
    for i=1:length(psths.V1_y)
    	cxy=psths.V1_y{i};
        thaly=psths.thal_y{i};
        cxy=cxy-nanmean(cxy(cxx>=params.spontWindow(1) & cxx<=params.spontWindow(2)));
        thaly=thaly-nanmean(thaly(cxx>=params.spontWindow(1) & cxx<=params.spontWindow(2)));
        psths.V1_y{i}=cxy;
        psths.thal_y{i}=thaly;
    end
end
if delayShift==1
    cxx=psths.V1_x{1};
    shiftInd=floor(0.003/(cxx(2)-cxx(1)));
    for i=1:length(psths.V1_y)
    	cxy=psths.V1_y{i};
        thaly=psths.thal_y{i};
        cxy=[cxy(shiftInd+1:end) zeros(1,shiftInd)];
        psths.V1_y{i}=cxy;
        psths.thal_y{i}=thaly;
    end
end

cxy=psths.V1_y{1};
cxx=psths.V1_x{1};
thaly=psths.thal_y{1};
% estimated_nsr=var(cxy(cxx>=params.spontWindow(1) & cxx<=params.spontWindow(2)))/5;
estimated_nsr=var(cxy(cxx>=params.spontWindow(1) & cxx<=params.spontWindow(2)))/5;
if test==1
%     taux=0:0.001:0.1;
%     tauy=-exp(-taux./0.01);
%     predCx=conv(tauy,thaly);
%     predCx=predCx(1:length(thaly));
%     tauguess=deconvwnr(predCx,thaly,estimated_nsr);
%     figure(); 
%     plot(0:0.001:0.001*(length(tauguess)-1),tauguess);
%     return

    taux=0:0.001:0.1;
    tauy=-exp(-taux./0.01);
    for i=1:length(psths.thal_y)
        thaly=psths.thal_y{i};
        cxy=conv(tauy,thaly);
        cxy=cxy(1:length(thaly));
        psths.V1_y{i}=cxy;
    end
end
if onlyStim==1
    cxy=cxy(cxx>=params.stimWindow(1) & cxx<=params.stimWindow(2));
    thaly=thaly(cxx>=params.stimWindow(1) & cxx<=params.stimWindow(2));
end
tauguess=deconvwnr(cxy,thaly,estimated_nsr);
estimate_taus=nan(length(params.useLEDcond),length(tauguess));
estimate_taus(1,:)=tauguess;

for i=2:length(params.useLEDcond)
    disp('calculating tau estimate');
	cxy=psths.V1_y{i};
    thaly=psths.thal_y{i};
    if onlyStim==1
        cxy=cxy(cxx>=params.stimWindow(1) & cxx<=params.stimWindow(2));
        thaly=thaly(cxx>=params.stimWindow(1) & cxx<=params.stimWindow(2));
    end
    tauguess=deconvwnr(cxy,thaly,estimated_nsr);
    estimate_taus(i,:)=tauguess;
end

% Align tau estimates -- because variable thalamus-cortex phase
% relationship, because sometimes you get On cells and sometimes you get
% Off cells
if doAlign==1 
    temp=downSampAv(estimate_taus(1,:),5);
    ds_est=nan(size(estimate_taus,1),length(temp));
    ds_est(1,:)=temp;
    for i=2:size(estimate_taus,1)
        ds_est(i,:)=downSampAv(estimate_taus(i,:),5);
    end
    aligned_est=nan(size(estimate_taus,1),2*size(estimate_taus,2));
    for i=1:size(ds_est,1)
%         [~,mind]=max(ds_est(i,350:450));
        [~,mind]=max(ds_est(i,395:405));
        startInd=size(estimate_taus,2)-(mind*5)-floor(size(estimate_taus,2)/2);
        aligned_est(i,startInd:startInd+size(estimate_taus,2)-1)=estimate_taus(i,:);
    end
    
    plotWindow=[3.5 4.5];
    tauguessx=0:0.001:0.001*(size(aligned_est,2)-1);
    if showFigs==1
        figure();
        o=0;
        usex=downSampAv(tauguessx(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),5);
        for i=1:size(aligned_est,1)
            curry=aligned_est(i,tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2));
            curry=downSampAv(curry,5);
            curry=curry-min(curry);
            plot(usex,curry+o);
            o=o+max(curry);
            hold on;
        end
    end
else
    tauguessx=0:0.001:0.001*(size(estimate_taus,2)-1);
    plotWindow=[1.5 2.5];
    aligned_est=estimate_taus;
end

tauguessy=nanmean(aligned_est,1);
tauguessse=nanstd(aligned_est,1)./sqrt(size(aligned_est,1));
if showFigs==1
    figure();
    plot(tauguessx(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),tauguessy(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)));
    hold on;
    plot(tauguessx(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),tauguessy(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2))+tauguessse(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),'Color',[0.5 0.5 0.5]);
    plot(tauguessx(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),tauguessy(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2))-tauguessse(tauguessx>=plotWindow(1) & tauguessx<=plotWindow(2)),'Color',[0.5 0.5 0.5]);
end

save([saveDirName '\estimate_taus.mat'],'aligned_est');
save([saveDirName '\tauguessx.mat'],'tauguessx');
save([saveDirName '\tauguessy.mat'],'tauguessy');
save([saveDirName '\tauguessse.mat'],'tauguessse');
end

function [varargout]=psth_wStdev(spikes,binsize,bsmooth,duration)
if nargin < 2
    binsize = 50; 
end

if nargin < 3
    bsmooth = 0;
end

if isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end

spiketimes = spikes.spiketimes;

binsize = binsize/1000;

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
centers = edges + diff(edges(1:2))/2;
if bsmooth
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end
varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),1);
end
