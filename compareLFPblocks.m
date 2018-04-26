function params=compareLFPblocks(expt)
global spikesDir physChannel photoChannel LFPdef_lowerCutoff LFPdef_upperCutoff BPLowerCutoff BPUpperCutoff dataDir;
% Obsolete version of wrapper for automatic experiment analysis

spikesDir='C:\Documents and Settings\Admin\My Documents\MATLAB\Data\Analyzed\SortedSpikes';
dataDir='C:\Documents and Settings\Admin\My Documents\MATLAB\Data\Rawdata';
physChannel=4;
photoChannel=2;
LFPdef_lowerCutoff=0; % in Hz   LFP definition
LFPdef_upperCutoff=200;  % in Hz   LFP definition
BPLowerCutoff=30; % in Hz   frequency band of interest
BPUpperCutoff=80; % in Hz   frequency band of interest

% Create output directory for figures
parentSaveDir=strcat('C:\Documents and Settings\Admin\My Documents\CompareLFP\',expt.name);
if ~exist(parentSaveDir,'dir')
    mkdir(parentSaveDir);
end

% Initialize parameters for analyses
params=getAnalysisParams(expt);

% For each block of constant stim. parameters, 
% check with user to be sure all the found DAQ files should be included, or
% get which should be excluded
analysisConfigs=[];
daqs1=[];
daqs2=[];
for i=1:length(params)
    paramsName=strcat('StimConfig',num2str(i));
    disp(strcat('Stim. Type:', params(i).stimType));
    disp(strcat('DAQ Files:', num2str(params(i).daqFiles)));
    exclude=input('Exclude the following DAQ files from this group. Enter answer in the format of an integer array (i.e., [1 2 3]).');
    useTheseDaqFiles=[];
    for j=params(i).daqFiles
        if any(j==exclude)
        else
            useTheseDaqFiles=[useTheseDaqFiles j];
        end
    end
    params(i).daqFiles=useTheseDaqFiles;
    % Analyze this stimulus parameter block separately
    % Passing in -100 to blockTrials1 or blockTrials2 causes function to
    % use all trials
    [savedIntoDir,analysisInfo]=analyzeLFPduringStim(expt,daqs1,daqs2,0,params(i),paramsName,[],[],parentSaveDir,1,-100,[]);
    analysisConfigs=[analysisConfigs; analysisInfo];
end

% Compare stimulus parameter blocks
stimConfig1=1;
paramName1=strcat('StimConfig',num2str(stimConfig1));
stimConfig2=2;
paramName2=strcat('StimConfig',num2str(stimConfig2));
[savedIntoDir,analysisInfo]=analyzeLFPduringStim(expt,daqs1,daqs2,1,params(stimConfig1),paramName1,params(stimConfig2),paramName2,1,parentSaveDir,-100,-100);

stimConfig1=1;
paramName1=strcat('StimConfig',num2str(stimConfig1));
stimConfig2=2;
paramName2=strcat('StimConfig',num2str(stimConfig2));
[savedIntoDir,analysisInfo]=analyzeLFPduringStim(expt,daqs1,daqs2,1,params(stimConfig1),paramName1,params(stimConfig2),paramName2,1,parentSaveDir,-100,-100);
    