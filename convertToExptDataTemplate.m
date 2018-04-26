function exptData=convertToExptDataTemplate(xpoints,y1,y2,params)

expt.name='KR_2012-03-30_M89A';
dataDir='Y:\New Acquisition Computer\';
saveDir='W:\Analysis Computer\Across Mice Finals\Anesth Evoked Silencing Across Mice\';

spontExpt=0;
checkLED=0;
% Get file names
if ~isempty(expt.name) && checkLED==1
    FileList = GetFileNames(expt,dataDir);
    expt.files.names = FileList;
else
    expt.files.names = [];
end

aliLED=0;
if isempty(params)
    params.useFileInd=[87:131];
    params.useStimcond=[1:8];
    params.useTrigger=[1:12];
    params.useLEDcond={[1]; [3]};
    params.trialDuration=6; % in s
    params.stimulusOn=[3 6]; % in s relative to trial onset
    params.stimulusType='DG';
    params.stimulusDuration=300; % in ms, stimulus duration before LED onset
    params.stimulusContrast=1;
    params.stimulusSize=1000; % in pixels, 1000 = full-field stimulus
    params.anesthOrAwake='anesth';
    params.anesthType='iso';
    params.LEDintensity=5;
end
params.LEDduration=1000; % in ms, duration of LED pulse -- maybe short to exclude photo-artifact
params.peakWait=0.05; % Time in s 'til response onset/peak
fracSuppParams.baseline=[0.6 0.999];
fracSuppParams.wait=0.04;

params.peakWindow=[params.stimulusOn(1)+params.peakWait params.stimulusOn(1)+(params.stimulusDuration/1000)-0.01];
fracSuppParams.ledWindow=[params.stimulusOn(1)+(params.stimulusDuration/1000)+fracSuppParams.wait params.stimulusOn(1)+(params.stimulusDuration/1000)+(params.LEDduration/1000)];

if spontExpt==0
    disp('DOING EVOKED');
else
    disp('DOING SPONT');
end

% Set up save directory for this mouse
saveDirName=[saveDir expt.name];
if exist(saveDirName,'dir')==0
    st=mkdir(saveDirName);
    if ~st
        disp('Could not create save directory');
    end
end

% Check LED alignment
if checkLED==1
    checkLEDonsets(expt,params.useFileInd,dataDir);
end

save([saveDirName '\paramsForSilencing.mat'],'params');
exptData.xpoints=xpoints;
exptData.ypoints1=y1;
exptData.ypoints2=y2;
exptData.numtrials1=nan;
exptData.numtrials2=nan;
save([saveDirName '\exptData.mat'],'exptData');

% Calculate fractional suppression for the summed unit PSTH and save
fracSupp=calcFractionalSupp_passInParams(exptData.xpoints,exptData.ypoints1,exptData.ypoints2,fracSuppParams,spontExpt);
save([saveDirName '\fracSupp_params.mat'],'fracSuppParams');
save([saveDirName '\fracSupp.mat'],'fracSupp');

% figure(); 
plot(exptData.xpoints,exptData.ypoints1,'Color','k');
hold on;
plot(exptData.xpoints,exptData.ypoints2,'Color','r');
Lonset=params.stimulusOn(1)+(params.stimulusDuration/1000);
line([Lonset Lonset],[min([exptData.ypoints1 exptData.ypoints2]) max([exptData.ypoints1 exptData.ypoints2])],'Color','b');
line([Lonset+(params.LEDduration/1000) Lonset+(params.LEDduration/1000)],[min([exptData.ypoints1 exptData.ypoints2]) max([exptData.ypoints1 exptData.ypoints2])],'Color','b');

end

% --- Subfunctions --- %
function FileList = GetFileNames(expt,DataDir)

% Get list of daq files for experiment in raw data directory
files = dir([DataDir expt.name '_' '*.daq']);
FileList = {files.name}';
if isempty(FileList)
    disp('No daq files found for this experiment name')
else
    % Need to reorder FileList (what's the best way to do this?)
    FilesFound = 0;
    i = 1;
    while(FilesFound < length(FileList))
        file = dir([DataDir expt.name '_' num2str(i) '.daq']);
        if ~isempty(file)
            FilesFound = FilesFound + 1;
            FileList(FilesFound) = {file.name}';
        end
        i = i+1;
        if i == 500
            disp('No daq files found for this experiment name');
        end
    end
    disp('The following daq files were found:')
    disp(FileList)
end
end