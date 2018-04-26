function analyzeResponseToStim(expt,params,currSaveDir,discardLastSweep,showPhoto,useLED,downSampFactor,nTrialsToUseForSpec,skipFirstTrial)
% Analyze the LFP and spiking responses to this stimulus presentation
% For a stimulus block, runs analysis, makes figures, saves figures, and
% saves photo-aligned data in .mat files
% 
% PARAMETERS:
% expt: the experiment structure currently being analyzed
% params: a structure with information about this stim. block (see
% getLFPAnalysisParams.m and allAnalyses.m)
% currSaveDir: directory into which to save figures and .mat data
% discardLastSweep: if 1, discards last sweep of each daq file
% showPhoto: if 1, show photodiode signal at the bottom of each graph
% useLED: if 1, separate stimulus conditions according to whether LED on or
% off
% downSampFactor: factor by which to down-sample data for speed
% nTrialsToUseForSpec: number of trials to average for spectrograms
% skipFirstTrial: if 1, skip the first trial of each daq file (need to fix
% this to actually skip first TWO trials of each daq file, because wrong
% LED pulse may occur on second trial)

global saveToDir2 spikesDir physChannel photoChannel LFPdef_lowerCutoff LFPdef_upperCutoff BPLowerCutoff BPUpperCutoff ledChannel photoAlignedDir LFPSaveDir;

% Get DAQ files to include in the analysis
daqFileNames=expt.files.names(params.daqFiles);
disp('Using these daq files');
disp(daqFileNames);

% Make output directory for units and LFP
saveToDir=strcat(currSaveDir,params.stimType,'');
if ~exist(saveToDir,'dir')
    mkdir(saveToDir);
end
saveToDir2=saveToDir;

% See if photodiode-aligned data already exists
% If so, load
% Else get data and align to photodiode
if ~exist(currSaveDir,'dir')
    mkdir(currSaveDir);
end
photoData=what(strcat(currSaveDir,'\PhotoAligned_Data'));
% Directory should contain photodiode-aligned LFPbySweep,
% bandPassedLFPbySweep, photodiodeBySweep and spikes
if isempty(photoData)
    alignData=1;
elseif isempty(photoData.mat)
    alignData=1;
elseif length(photoData.mat)<4
    alignData=1;
else
    alignData=0;
end
if alignData==1   
    % Get LFP -- LFP will be same for every unit
    [LFPbySweep,LFP_Fs,bandPassedLFPbySweep,photodiodeBySweep,ledAv,ledBySweep,gotDaqs,completeSweeps]=getLFPbySweep(daqFileNames,params.samplingRate,LFPdef_lowerCutoff,LFPdef_upperCutoff,BPLowerCutoff,BPUpperCutoff,downSampFactor,physChannel,photoChannel,ledChannel);
    if any(size(LFPbySweep)==0)
        disp('No data in these daq files.');
    end
    daqFileNames=daqFileNames(gotDaqs);
    params.daqFiles=params.daqFiles(gotDaqs); 
    
    % Get stimulus condition corresponding to each trial in LFPbySweep
    stimsForSweeps=[];
    for i=params.daqFiles 
        stimsForSweeps=[stimsForSweeps expt.sweeps.stimcond(expt.sweeps.fileInd==i)];
    end
    
    % KR hack
    %stimsForSweeps=ones(1,length(stimsForSweeps));
    
    
    % Either rely on LED condition file to have the correct LED condition
    % for each sweep (which it sometimes will not, particularly if
    % acquiring short sweeps, or if the LED is ON for every trial)
    % OR read in the LED condition from the appropriate analog input
    % channel
    % Get LED condition from analog input channel
    [ledForSweeps,mForLED]=getLEDConditions(ledBySweep,LFP_Fs);
    if mForLED<0.0002
        % Use LED condition file
        ledForSweeps=[];
        for i=params.daqFiles
            ledForSweeps=[ledForSweeps expt.sweeps.led(expt.sweeps.fileInd==i)];
        end
    end   
    %clear ledBySweep
   
    
    % If no LED data, assume LED is off for all trials
    if isempty(ledForSweeps)
        ledForSweeps=zeros(1,size(LFPbySweep,1));
    end
    % If no stim. conditions, assume all trials have same stim. condition
    if isempty(stimsForSweeps)
        stimsForSweeps=ones(1,size(LFPbySweep,1));
    end
    if ~any(~(isnan(ledForSweeps))) % All NaN
        ledForSweeps=zeros(1,size(LFPbySweep,1));
    end
    if ~any(~(isnan(stimsForSweeps))) % All NaN
        stimsForSweeps=ones(1,size(LFPbySweep,1));
    end
    
    if ~useLED % KR hack
        ledForSweeps=zeros(1,length(stimsForSweeps));
    end
    
    % KR hack
    % If using artifact subtraction, subtract now from trials with LED on
%     a=load('E:\Results\GammaLFPs\KR_2010-08-23\LED Artifacts\AI76_LEDartifact.mat');
%     artifact=a.artifact;
%     LEDonset=0.8;
%     %LFPbySweep(:,floor(0.5*LFP_Fs)+1:floor(0.5*LFP_Fs)+1+2560-1)=LFPbySweep(:,floor(0.5*LFP_Fs)+1:floor(0.5*LFP_Fs)+1+2560-1)-(ones(size(LFPbySweep,1),1)*artifact);
%     LFPbySweep(ledForSweeps>0,floor(LEDonset*LFP_Fs)+1:floor(LEDonset*LFP_Fs)+1+length(artifact)-1)=LFPbySweep(ledForSweeps>0,floor(LEDonset*LFP_Fs)+1:floor(LEDonset*LFP_Fs)+1+length(artifact)-1)-(ones(length(find(ledForSweeps>0)),1)*artifact);
%     
    % Unit specific-analyses for all units in this experiment
    try
        APs=load(strcat(spikesDir,'\',expt.sort.trode.spikesfile));
        spikes=APs.spikes;
        spikes=filtspikes(spikes,1,'fileInd',params.daqFiles);
        units=unique(spikes.assigns);
    catch
        disp('Could not load spikes for this experiment.');
        spikes.trials=[];
        units=[];
    end
    
    % Align to photodiode
    if ~any(size(LFPbySweep)==0)
        [spikes,LFPbySweep,LFP_Fs,bandPassedLFPbySweep,photodiodeBySweep,ledBySweep,params]=photoAlignData(spikes,LFPbySweep,LFP_Fs,bandPassedLFPbySweep,photodiodeBySweep,ledBySweep,params);
    end
        
    % Save data aligned to photodiode
    photoAlignedDir=strcat(currSaveDir,'\PhotoAligned_Data');
    if ~exist(photoAlignedDir,'dir')
        mkdir(photoAlignedDir);
    end
    save(strcat(photoAlignedDir,'\spikes.mat'),'spikes');
    save(strcat(photoAlignedDir,'\LFPbySweep.mat'),'LFPbySweep','stimsForSweeps');
    save(strcat(photoAlignedDir,'\bandPassedLFPbySweep.mat'),'bandPassedLFPbySweep');
    save(strcat(photoAlignedDir,'\photodiodeBySweep.mat'),'photodiodeBySweep');
    save(strcat(photoAlignedDir,'\LFP_Fs.mat'),'LFP_Fs');
    save(strcat(photoAlignedDir,'\stimsForSweeps.mat'),'stimsForSweeps');
    save(strcat(photoAlignedDir,'\ledForSweeps.mat'),'ledForSweeps');
    save(strcat(photoAlignedDir,'\ledAv.mat'),'ledAv');
    save(strcat(photoAlignedDir,'\ledBySweep.mat'),'ledBySweep');
    p=params.ONstart;
    save(strcat(photoAlignedDir,'\ONstart.mat'),'p');
    save(strcat(photoAlignedDir,'\completeSweeps.mat'),'completeSweeps');
else
    % Just read in mat files
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\spikes.mat'));
    spikes=a.spikes;
    %spikes=filtspikes(spikes,1,'fileInd',params.daqFiles); % KR hack! for
    %LFP analysis
    if isfield(spikes,'assigns')
        units=unique(spikes.assigns);
    else
        disp('No units loaded for this experiment.');
        units=[];
    end
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\LFPbySweep.mat'));
    LFPbySweep=a.LFPbySweep;
    stimsForSweeps=a.stimsForSweeps;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\bandPassedLFPbySweep.mat'));
    bandPassedLFPbySweep=a.bandPassedLFPbySweep;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\photodiodeBySweep.mat'));
    photodiodeBySweep=a.photodiodeBySweep;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\LFP_Fs.mat'));
    LFP_Fs=a.LFP_Fs;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\stimsForSweeps.mat'));
    stimsForSweeps=a.stimsForSweeps;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\ledForSweeps.mat'));
    ledForSweeps=a.ledForSweeps;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\ledAv.mat'));
    ledAv=a.ledAv;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\ONstart.mat'));
    params.ONstart=a.p;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\completeSweeps.mat'));
    completeSweeps=a.completeSweeps;
end  

if any(size(LFPbySweep)==0)
    return % No data for the set of daq files
end

% if exist('E:\Results\GammaLFPs\KR_2010-08-02\currLEDForTrials.mat','file')
%     a=load('E:\Results\GammaLFPs\KR_2010-08-02\currLEDForTrials.mat');
%     ledForSweeps=a.ledForSweeps;
% end


% Make spectrogram - hacked in
% sig=cell(150,1);
% these=randperm(size(LFPbySweep,1));
% for i=1:150
%     sig{i}=LFPbySweep(these(i),floor((params.ONstart+params.ONlength-0.35)*LFP_Fs):end)';
% end
% trySpecgram(sig,LFP_Fs);
% 'hi'

% Iterate through units
% Define a structure for saving info. on unit-specific analyses
unitAnalysisInfo=struct();
unitAnalysisInfo.definedRedTrials=[];
unitAnalysisInfo.definedBlueTrials=[];
unitAnalysisInfo.redTrial_spikingThresh=[];
unitAnalysisInfo.responseQuantMethod=[];
unitAnalysisInfo.stimConds_colorCode=[];
unitAnalysisInfo.comparedRedTrials=[];
unitAnalysisInfo.comparedBlueTrials=[];
%for u=units
for u=[]
    unitSaveDir=strcat(saveToDir,'\','Assigns',num2str(u));
    if ~exist(unitSaveDir,'dir')
        mkdir(unitSaveDir);
    end
    % Analyze this unit
    [unitFigs,unitAnalysisInfo]=analyzeUnit(expt,filtspikes(spikes,0,'assigns',u),params,LFPbySweep,bandPassedLFPbySweep,LFP_Fs,unitAnalysisInfo);
    % LFP Analyses
    [LFPfigs,usedRed,usedBlue]=analyzeLFPSweeps(LFPbySweep,bandPassedLFPbySweep,LFP_Fs,unitAnalysisInfo.comparedRedTrials,unitAnalysisInfo.comparedBlueTrials,params); 
    
    % Save figures from this analysis and the analysis details used to
    % produce them
    saveas(unitFigs(1),strcat(unitSaveDir,'\','PSTH.bmp'));
    saveas(unitFigs(2),strcat(unitSaveDir,'\','PSTH_separateStimSpecific.bmp'));
    saveas(unitFigs(3),strcat(unitSaveDir,'\','polarPlots.bmp'));
    saveas(unitFigs(4),strcat(unitSaveDir,'\','rasters.bmp'));
    saveas(unitFigs(5),strcat(unitSaveDir,'\','spikeTriggeredLFP.bmp'));
    saveas(unitFigs(6),strcat(unitSaveDir,'\','spikeTriggeredLFP_diffOrientations.bmp'));
    saveas(LFPfigs(1),strcat(unitSaveDir,'\','stimulusTriggeredLFP.bmp')); 
    
    saveas(unitFigs(1),strcat(unitSaveDir,'\','PSTH.fig'));
    saveas(unitFigs(2),strcat(unitSaveDir,'\','PSTH_separateStimSpecific.fig'));
    saveas(unitFigs(3),strcat(unitSaveDir,'\','polarPlots.fig'));
    saveas(unitFigs(4),strcat(unitSaveDir,'\','rasters.fig'));
    saveas(unitFigs(5),strcat(unitSaveDir,'\','spikeTriggeredLFP.fig'));
    saveas(unitFigs(6),strcat(unitSaveDir,'\','spikeTriggeredLFP_diffOrientations.fig'));
    saveas(LFPfigs(1),strcat(unitSaveDir,'\','stimulusTriggeredLFP.fig'));  
    save(strcat(unitSaveDir,'\','AnalysisDetails.mat'),'unitAnalysisInfo');
    LFPanalysisInfo.usedRed=usedRed;
    LFPanalysisInfo.usedBlue=usedBlue;
    save(strcat(unitSaveDir,'\','LFPAnalysisDetails.mat'),'LFPanalysisInfo');
end  

% Compare LFP responses by stimulus condition
LFPSaveDir=strcat(saveToDir,'\LFP');
if ~exist(LFPSaveDir,'dir')
    mkdir(LFPSaveDir);
end

stimsForSweeps=stimsForSweeps';
ledForSweeps=ledForSweeps';

% Discard truncatedSweeps
if 1
    LFPbySweep=LFPbySweep(logical(completeSweeps),:);
    stimsForSweeps=stimsForSweeps(logical(completeSweeps),:);
    ledForSweeps=ledForSweeps(logical(completeSweeps),:);
    bandPassedLFPbySweep=bandPassedLFPbySweep(logical(completeSweeps),:);
end

% If discardLastSweep==1, discard last sweep
if discardLastSweep==1 
    LFPbySweep=LFPbySweep(1:end-1,:);
    stimsForSweeps=stimsForSweeps(1:end-1,:);
    ledForSweeps=ledForSweeps(1:end-1,:);
    bandPassedLFPbySweep=bandPassedLFPbySweep(1:end-1,:);
end
% If skipFirstTrial==1, discard first sweep
if skipFirstTrial==1 
    LFPbySweep=LFPbySweep(2:end,:);
    stimsForSweeps=stimsForSweeps(2:end,:);
    ledForSweeps=ledForSweeps(2:end,:);
    bandPassedLFPbySweep=bandPassedLFPbySweep(2:end,:);
end

if any(size(LFPbySweep)==0)
    disp('Only 1 sweep in this set of daq files. May be truncated. Not using.');
    return % No data for the set of daq files
end

% Get rid of trials with stimulus or led condition == nan (not useful)
bothInds=~(isnan(stimsForSweeps)|isnan(ledForSweeps));
stimsForSweeps=stimsForSweeps(bothInds);
ledForSweeps=ledForSweeps(bothInds);
LFPbySweep=LFPbySweep(bothInds,:);
bandPassedLFPbySweep=bandPassedLFPbySweep(bothInds,:);
% sIndex=stimsForSweeps(bothInds);
% lIndex=ledForSweeps(bothInds);
% sIndex=isnan(stimsForSweeps);
% lIndex=isnan(ledForSweeps);
% if length(sIndex)~=length(lIndex) % In case no stimuli
%     disp('Number of stimulus and LED conditions do not match.');
%     sIndex=ones(length(lIndex),1);
%     params.Var1_values=1;
%     params.Var2_values=1;
% end
% LFPbySweep=LFPbySweep(~(sIndex | lIndex),:);
% stimsForSweeps=stimsForSweeps(~(sIndex | lIndex),:);
% ledForSweeps=ledForSweeps(~(sIndex | lIndex),:);
% bandPassedLFPbySweep=bandPassedLFPbySweep(~(sIndex | lIndex),:);

photoSignal=mean(photodiodeBySweep(1:end-1,:),1);
clear photodiodeBySweep

if any(size(LFPbySweep)==0)
    disp('No data!');
    return
end

% LFP filtered 0 to 300 Hz
[lineColors,colorMapCode]=makeColorCode(params);
[LFPfigs,analysisInfo]=analyzeLFPbyStimCond(LFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'average',[],[],colorMapCode,0,50,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
for i=1:length(LFPfigs)
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_0to300Hz',num2str(i),'.bmp'));
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_0to300Hz',num2str(i),'.fig'));
end
%saveas(LFPfigs(1),strcat(LFPSaveDir,'\',''));


% Get alpha and beta filtered LFP data
% alphaLFPbySweep=bandPassLFP(LFPbySweep,LFP_Fs,8,12,1);
% betaLFPbySweep=bandPassLFP(LFPbySweep,LFP_Fs,12,30,1);

clear LFPbySweep

% LFP filtered 30 to 80 Hz
[LFPfigs,analysisInfo]=analyzeLFPbyStimCond(bandPassedLFPbySweep,stimsForSweeps,ledForSweeps,'initial',LFP_Fs,params,[1 0 0],'average',[],[],colorMapCode,1,50,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
for i=1:length(LFPfigs)
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80Hz',num2str(i),'.bmp'));
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80Hz',num2str(i),'.fig'));
end

% Gamma-filtered LFP integral
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);
[LFPfigs,analysisInfo]=analyzeLFPbyStimCond(bandPassedLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,50,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
for i=1:length(LFPfigs)
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80HzIntegral',num2str(i),'.bmp'));
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80HzIntegral',num2str(i),'.fig'));
end
% Smoothed gamma-filtered LFP integral
for i=1:size(bandPassedLFPbySweep,1)
    bandPassedLFPbySweep(i,:)=smooth(bandPassedLFPbySweep(i,:)',250)'; % Approx. 9 ms window for smoothing, given Fs=32000
end
[LFPfigs,analysisInfo]=analyzeLFPbyStimCond(bandPassedLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,100,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
for i=1:length(LFPfigs)
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80HzSmoothedIntegral',num2str(i),'.bmp'));
    saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_30to80HzSmoothedIntegral',num2str(i),'.fig'));
end
clear bandPassedLFPbySweep
return
% % LFP filtered 8 to 12 Hz
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(alphaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'average',[],[],colorMapCode,1,100,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12Hz',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12Hz',num2str(i),'.fig'));
% end
% % % Alpha-filtered LFP integral
% alphaLFPbySweep(alphaLFPbySweep<0)=-alphaLFPbySweep(alphaLFPbySweep<0);
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(alphaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,100,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12HzIntegral',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12HzIntegral',num2str(i),'.fig'));
% end
% % % Smoothed alpha-filtered LFP integral
% for i=1:size(alphaLFPbySweep,1)
%     alphaLFPbySweep(i,:)=smooth(alphaLFPbySweep(i,:)',1000)'; % Approx. 1/2 alpha period, given Fs=32000
% end
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(alphaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,500,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12HzSmoothedIntegral',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_8to12HzSmoothedIntegral',num2str(i),'.fig'));
% end
% clear alphaLFPbySweep
% 
% % % LFP filtered 12 to 30 Hz
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(betaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'average',[],[],colorMapCode,1,50,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30Hz',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30Hz',num2str(i),'.fig'));
% end
% % Beta-filtered LFP integral
% betaLFPbySweep(betaLFPbySweep<0)=-betaLFPbySweep(betaLFPbySweep<0);
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(betaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,100,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30HzIntegral',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30HzIntegral',num2str(i),'.fig'));
% end
% % Smoothed beta-filtered LFP integral
% for i=1:size(betaLFPbySweep,1)
%     betaLFPbySweep(i,:)=smooth(betaLFPbySweep(i,:)',500)'; 
% end
% [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(betaLFPbySweep,stimsForSweeps,ledForSweeps,'mean',LFP_Fs,params,[1 0 0],'add',[],[],colorMapCode,1,250,showPhoto,useLED,photoSignal,ledAv,0,nTrialsToUseForSpec);
% for i=1:length(LFPfigs)
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30HzSmoothedIntegral',num2str(i),'.bmp'));
%     saveas(LFPfigs(i),strcat(LFPSaveDir,'\','LFPFigure_12to30HzSmoothedIntegral',num2str(i),'.fig'));
% end