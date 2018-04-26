function [ONtrials,OFFtrials,useAssigns,wwidths,cellPhases,isOFF,psths,unitAmps,freqsMatched,sameFreqSet,togAvResponseStim,unitNotFollowingAmps]=separateONandOFFcells(spikes,useAssigns,ONtrials,OFFtrials)

bin=1; % in ms
fileInd=[32:79];
trialDuration=10.5; % in s
stimWindow=[0.5 9.5];
dLGNdiscard=0; % will discard thin units
testFS=0; % ask user whether unit is FS
Fs=25000;
plotMedian=0;
% params.tapers=[15 18];
params.tapers=[10 18];
params.pad=0;
params.fpass=[0.5 70];
params.Fs=1000/bin;
params.trialave=1;   
downSampForFig=10;
freqbinForXcorr=3; % in Hz

if isempty(ONtrials)
    if dLGNdiscard==1
        newUseAssigns=[];
        wwidths=zeros(1,length(useAssigns));
        for i=1:length(useAssigns)
            use_a_inds=spikes.assigns==useAssigns(i);
            theseEvChs=spikes.info.detect.event_channel(use_a_inds);
            topEvCh=mode(theseEvChs);
            wvfrmsEvCh=spikes.waveforms(:,:,topEvCh);
            wvfrm=mean(wvfrmsEvCh(use_a_inds,:),1);
            wwidths(i)=classifyUnitWaveformWidth(wvfrm,3.2*10^-4,Fs);
            if wwidths(i)>=0.22*10^-3
                newUseAssigns=[newUseAssigns useAssigns(i)];
            end
        end
        useAssigns=newUseAssigns;
    elseif testFS==1
        newUseAssigns=[];
        wwidths=zeros(1,length(useAssigns));
        f=figure(); 
        movegui(f,'south'); 
        for i=1:length(useAssigns)
            use_a_inds=spikes.assigns==useAssigns(i);
            theseEvChs=spikes.info.detect.event_channel(use_a_inds);
            topEvCh=mode(theseEvChs);
            wvfrmsEvCh=spikes.waveforms(:,:,topEvCh);
            wvfrm=mean(wvfrmsEvCh(use_a_inds,:),1);
            plot(wvfrm);
            isFSans=questdlg('Is unit fast-spiking (FS)?');
            if strcmp(isFSans,'No')
                FS=0;
            elseif strcmp(isFSans,'Yes')
                FS=1;
            elseif strcmp(isFSans,'Cancel')
                return
            end
            if FS==0
                newUseAssigns=[newUseAssigns useAssigns(i)];
            end
        end
        useAssigns=newUseAssigns;
    else
        wwidths=[];
    end
    
    % Input stimulus
    stim_x=linspace(0,3,3*(1000/bin));
    stim_y_temp=chirp(stim_x,1,3,30,'logarithmic')+0.5;
    stim_y=[stim_y_temp fliplr(stim_y_temp) stim_y_temp];
    
    % Get unit psths
    [x,psths]=getTrialByTrialUnitPSTH(spikes,useAssigns,trialDuration,bin,fileInd);
    
    % For each unit, choose stimulus frequency that drives largest amplitude
    % response in unit, then see whether, in that cycle at that frequency, more
    % spikes occur in OFF phase or ON phase
    % Highest stimulus frequency is 30 Hz, so bin at 10 ms to get amplitude
    % forAmpStim=downSampAv(stim_y,floor(10/bin));
    % [~,stimPeaks]=findpeaks(forAmpStim);
    
    takeFraction=0.95;
    nTrials=10;
    [Sstim,tstim,fstim]=mtspecgrampb(stim_y',[1.5 0.05],params);
    chirpfreq=1*((30/1)^(1/3)).^stim_x;
    chirpfreq=[chirpfreq fliplr(chirpfreq) chirpfreq];
    stim_x=linspace(0,9,9*(1000/bin));
    useTheseInds=zeros(size(Sstim,1),size(Sstim,2));
    matchedChirpfreq=zeros(1,length(tstim));
    for i=1:length(tstim)
        currchirpfreq=mean([chirpfreq(find(stim_x>=tstim(i),1,'first')) chirpfreq(find(stim_x<=tstim(i),1,'last'))]);
        matchedChirpfreq(i)=currchirpfreq;
        fstimInds=[find(fstim<=currchirpfreq,1,'last') find(fstim>=currchirpfreq,1,'first')];
        useTheseInds(i,fstimInds)=[1 1];
    end
    doneMatched=zeros(size(matchedChirpfreq));
    freqsMatched=[];
    for i=1:length(matchedChirpfreq)
        if doneMatched(i)==1
            break
        end
        currmatched=find(matchedChirpfreq<matchedChirpfreq(i)+0.1 & matchedChirpfreq>matchedChirpfreq(i)-0.1);
        sameFreqSet{i}=currmatched;
        freqsMatched(i)=matchedChirpfreq(i);
        doneMatched(currmatched)=1;
    end
    isOFF=zeros(1,length(psths));
    cellPhases=zeros(1,length(psths));
    followingStim=Sstim.*useTheseInds;
    avResponseStim=sum(followingStim,2);
    togAvResponseStim=nan(1,length(sameFreqSet));
    for j=1:length(sameFreqSet)
        togAvResponseStim(j)=nanmean(avResponseStim(sameFreqSet{j}));
    end
%     ONwindows=[397 670; 879 969; 5124 5331; 5610 6002; 6398 6667; 6877 6969];
%     OFFwindows=[1 397; 670 879; 5032 5124; 5331 5610; 6002 6398; 6667 6877];
    [~,f]=findpeaks(stim_y);
    OFFstarts=[1 f];
    [~,f]=findpeaks(-stim_y);
    ONstarts=[f];
    OFFwindows=[OFFstarts' ONstarts'];
    ONwindows=[ONstarts(1:end-1)' OFFstarts(2:end)'];
    unitAmps=zeros(length(psths),length(sameFreqSet));
    for i=1:length(psths)
        disp(i);
        currcell=psths{i};
        currcell=currcell(:,x>=stimWindow(1) & x<=stimWindow(2));
        %     [S,t,f]=mtspecgrampb(currcell',[1 0.1],params);
        [S,t,f]=mtspecgrampb(nanmean(currcell,1)',[1.5 0.05],params);
        followingS=S.*useTheseInds;
        avResponse=sum(followingS,2);
        togAvResponse=nan(1,length(sameFreqSet));
        notUseTheseInds=useTheseInds;
        notUseTheseInds(notUseTheseInds==0)=2;
        notUseTheseInds(notUseTheseInds==1)=0;
        notUseTheseInds(notUseTheseInds==2)=1;
        otherS=S.*notUseTheseInds;
        notAvResponse=nanmean(otherS,2);
        notTogAvResponse=nan(1,length(sameFreqSet));
        for j=1:length(sameFreqSet)
            togAvResponse(j)=nanmean(avResponse(sameFreqSet{j}));
            notTogAvResponse(j)=nanmean(notAvResponse(sameFreqSet{j}));
        end
%         unitAmps(i,:)=sqrt(togAvResponse./togAvResponseStim);
        unitAmps(i,:)=sqrt(togAvResponse);
        unitNotFollowingAmps(i,:)=sqrt(notTogAvResponse);
        % Find frequency that drives the maximal response in this cell
        [~,maxInd]=max(togAvResponse(freqsMatched<=30)./togAvResponseStim(freqsMatched<=30));
        bestFreq=freqsMatched(maxInd);
        bestFreq=1.5;
        % Find the phase of this cell's response at this best frequency --
        % from cross-correlation
        timeIndsToUse=find(chirpfreq>=bestFreq-(freqbinForXcorr/2) & chirpfreq<=bestFreq+(freqbinForXcorr/2));
        c=nanmean(currcell,1);
%         zeroc=zeros(size(c));
%         zeroc(timeIndsToUse)=c(timeIndsToUse)-mean(c(timeIndsToUse));
%         zerostim=zeros(size(stim_y));
%         zerostim(timeIndsToUse)=stim_y(timeIndsToUse)-mean(stim_y(timeIndsToUse));
        zeroc=c;
        zerostim=stim_y;
        ONcount=0;
        ONinds=0;
        for j=1:size(ONwindows,1)
            ONcount=ONcount+sum(zeroc(ONwindows(j,1):ONwindows(j,2)));
            ONinds=ONinds+length(ONwindows(j,1):ONwindows(j,2));
        end
        OFFcount=0;
        OFFinds=0;
        for j=1:size(OFFwindows,1)
            OFFcount=OFFcount+sum(zeroc(OFFwindows(j,1):OFFwindows(j,2)));
            OFFinds=OFFinds+length(OFFwindows(j,1):OFFwindows(j,2));
        end
        cellPhases(i)=ONcount/ONinds;
%         if ONcount/ONinds>OFFcount/OFFinds
        if ONcount/ONinds>3.5
            isOFF(i)=0;
        else
            isOFF(i)=1;
        end
        
%         [cc,cc_lags]=xcorr(zeroc,zerostim,[],'none');
%         bestPeriod=1/bestFreq;
%         cc=cc(cc_lags>=-bestPeriod*(1000/bin) & cc_lags<=bestPeriod*(1000/bin));
%         cc_lags=cc_lags(cc_lags>=-bestPeriod*(1000/bin) & cc_lags<=bestPeriod*(1000/bin));
%         [~,maxccInd]=max(cc);
%         useLag=cc_lags(maxccInd);
%         periodInLagInds=bestPeriod*(1000/bin);
%         cellPhases(i)=useLag/periodInLagInds;
%         if cellPhases(i)>=-0.6 & cellPhases(i)<0.2
%             isOFF(i)=0;
%         elseif cellPhases(i)>=0.2 & cellPhases(i)<=1
%             isOFF(i)=1;
%         end

%         if useLag>=0 && useLag<periodInLagInds/2
%             isOFF(i)=0;
%         elseif useLag>=periodInLagInds/2 && useLag<=periodInLagInds
%             isOFF(i)=1;
%         elseif useLag>=-periodInLagInds/2 && useLag<0
%             isOFF(i)=1;
%         elseif useLag>=-periodInLagInds && useLag<-periodInLagInds/2
%             isOFF(i)=0;
%         end
        
        % Find the phase of this cell's response at this best frequency --
        % from COHERENCE PHI
%         bootcell=zeros(nTrials,length(nanmean(currcell,1)));
%         for j=1:nTrials
%             r=randi([1 size(currcell,1)],1,floor(size(currcell,1)*takeFraction));
%             bootcell(j,:)=nanmean(currcell(r,:),1);
%         end
%         [C,phi,~,~,~,f_coh]=coherencypb(repmat(stim_y',1,nTrials),bootcell',params,0);
%         bestPhase=nanmean([phi(find(f_coh<=bestFreq,1,'last')) phi(find(f_coh>=bestFreq,1,'first'))]);
%         cellPhases(i)=bestPhase;
%         if bestPhase>=-pi/2 & bestPhase<=pi/2
%             isOFF(i)=0;
%         else
%             isOFF(i)=1;
%         end
        
        
        %     forAmp=nanmean(downSampMatrix(currcell,floor(10/bin)),1);
        %     [~,maxInd]=max(forAmp);
        %     startCycle=find(stimPeaks<=maxInd,1,'last');
        %     endCycle=find(stimPeaks>maxInd,1,'first');
        %     midCycle=floor(startCycle+(endCycle-startCycle)/2);
        %     OFFphaseSpikes=sum(forAmp(startCycle:midCycle));
        %     ONphaseSpikes=sum(forAmp(midCycle+1:endCycle));
        %     if OFFphaseSpikes>ONphaseSpikes
        %         isOFF(i)=1;
        %     end
    end
    
    % Get trial-by-trial ON cell and OFF cell averages
    ONtrials=zeros(size(psths{1}));
    OFFtrials=zeros(size(psths{1}));
    for i=1:length(psths)
        currcell=psths{i};
        if isOFF(i)==1
            OFFtrials=OFFtrials+currcell;
        else
            ONtrials=ONtrials+currcell;
        end
    end
    OFFtrials=OFFtrials/sum(isOFF);
    ONtrials=ONtrials/sum(~isOFF);
else
    x=linspace(0,10.5,size(ONtrials,2));
end

% Plot ON and OFF averages
figure(); 
if plotMedian==1
    plot(downSampAv(x,downSampForFig),downSampAv(median(ONtrials,1),downSampForFig),'Color','r');
else
    plot(downSampAv(x,downSampForFig),downSampAv(nanmean(ONtrials,1),downSampForFig),'Color','r');
end
hold on;
if plotMedian==1
    plot(downSampAv(x,downSampForFig),downSampAv(median(OFFtrials,1),downSampForFig),'Color','b');
else
    plot(downSampAv(x,downSampForFig),downSampAv(nanmean(OFFtrials,1),downSampForFig),'Color','b');
end

end

function [x,psths_t]=getTrialByTrialUnitPSTH(spikes,allAssigns,trialDuration,bin,fileInd)

useStimcond={[1:1000]};
useLED=[];

% Get trials for each unit
unitByUnitTrials=cell(1,length(allAssigns));
unitByUnitStimcond=cell(1,length(allAssigns));
unitByUnitLED=cell(1,length(allAssigns));
unitByUnitFileInd=cell(1,length(allAssigns));
for i=1:length(allAssigns)
    unitByUnitTrials{i}=unique(spikes.sweeps.trials);
    unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
    unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
    unitByUnitFileInd{i}=spikes.sweeps.fileInd(unique(spikes.sweeps.trials));
end

psths_t=cell(length(allAssigns),length(useStimcond));
x=[];
for i=1:length(allAssigns)
    disp(i);
    for j=1:length(useStimcond)
        useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        currFileIndCon=unitByUnitFileInd{i};
        if ~isempty(useLED) && ~isempty(fileInd)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED) & ismember(currFileIndCon,fileInd));
        elseif ~isempty(useLED) && isempty(fileInd)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED));
        elseif ~isempty(fileInd)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currFileIndCon,fileInd));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        [~,~,~,x,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
    end
end

end

function [varargout]=psth_wStd_trialByTrial(spikes,binsize,bsmooth,duration,nTrials,theseTrials)

if nargin < 2
    binsize = 50; 
end
if nargin < 3
    bsmooth = 1;
end
% Set duration and number of trials
if ~isempty(nTrials)
    numtrials=nTrials;
elseif isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

nsForStdev=zeros(numtrials,size(n,2));
if ~isempty(theseTrials)
    allTrials=theseTrials;
else
    allTrials=unique(spikes.trials);
end
if length(allTrials)~=numtrials
    if ~isempty(theseTrials)
        allTrials=theseTrials;
    elseif length(spikes.sweeps.trials)==numtrials
        allTrials=spikes.sweeps.trials;
    else
        disp('Needed to fill in trials -- be sure you are using contiguous daq files');
        allTrials=min(unique(spikes.trials)):max(unique(spikes.trials));
    end
end      
for i=1:length(allTrials)
    cspikes=filtspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
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
% varargout{6} = 1*std(nsForStdev(:,1:end-1),0,1);
varargout{6} = nsForStdev(:,1:end-1);
end

function [varargout]=calcMeanAndStdForUnit(spikes,window,theseTrials)

binsize=1; % in ms
if isempty(theseTrials)
    disp('theseTrials should not be empty in calcMeanAndStdForUnit');
end
numtrials=length(theseTrials);
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges=window(1):binsize:window(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

nsForStdev=zeros(numtrials,size(n,2));
allTrials=theseTrials;
for i=1:length(allTrials)
    cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
s=std(mean(nsForStdev,2),0,1);
if all(isnan(n))
    n = 0;
end
varargout{1} = mean(mean(nsForStdev,2));
varargout{2} = s;
varargout{3} = mean(nsForStdev,2);
end

