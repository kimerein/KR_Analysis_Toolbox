function allAnalyses(expt,stimBlocks,excludeTheseDaqs)
% Wrapper function for automated analysis of an experiment
% Pass in 
% 1. the desired experiment structure, expt
% 2. stimBlocks, which is a matrix specifying which 
% daq files (as referenced by their indices in expt.files.names, this index 
% is called the "analysis index" for each daq file)
% should be included in each separate stimulus blocks
% Each stimulus block is analyzed separately and results (including figures
% and organized + photo-aligned data) are saved to the directory for that 
% experiment name. If stimBlocks is [], stimulus blocks are inferred from
% data in the expt structure.
% First column of stimBlocks should give the index of the first daq file to
% include in each stim. blocks; second column should give the index of the
% last daq file to include in each stim. block. 
% The number of rows of stimBlocks should be the number of separate stim.
% blocks for analysis. 
% e.g., 
% stimBlocks = [1 2; 3 7; 8 10] makes three separate blocks for analysis
% with daq files #s (analysis indices) 1+2 in the first block, 
% 3,4,5,6,7 in the second block
% and 8,9,10 in the third block
% 3. excludeTheseDaqs specifies which daq files (according to analysis
% indices) to exclude from each stimBlocks
% excludeTheseDaqs is a cell array with each cell containing a vector of 
% daq file analysis indices. The daq files in excludeTheseDaqs will be
% removed from the analysis. 
% e.g., 
% excludeTheseDaqs = {[2]; [5 6]; [9]} plus stimBlocks = [1 2; 3 7; 8 10]
% will give three analysis blocks containing, respectively, data from
% daq files (by analysis index) 
% 1
% 3 4 7
% 8 10

global spikesDir ledChannel physChannel photoChannel LFPdef_lowerCutoff LFPdef_upperCutoff BPLowerCutoff BPUpperCutoff dataDir;
global waveletDir

spikesDir='E:\MATLAB\Data\Analyzed\SortedSpikes'; % Directory containing sorted spikes
%spikesDir='C:\Documents and Settings\Admin\My Documents\MATLAB\Data\Analyzed\SortedSpikes';
dataDir='E:\MATLAB\Data\RawData'; % Directory containing raw .daq data
%dataDir='C:\Documents and Settings\Admin\My Documents\MATLAB\Data\Rawdata';
physChannel=11; % Matlab index for the physiology channel
photoChannel=6; % Matlab index for the photodiode channel
ledChannel=4; % Matlab index for the LED channel

% physChannel=6;
% photoChannel=8;
% ledChannel=9;

LFPdef_lowerCutoff=0.0001; % in Hz   LFP definition
LFPdef_upperCutoff=200;  % in Hz   LFP definition
BPLowerCutoff=30; % in Hz   frequency band of interest
BPUpperCutoff=80; % in Hz   frequency band of interest
discardLastSweep=1; % Discard last sweep of each daq file, if 1
discardFirstSweep=1; % Discared first sweep of each daq file, if 1 -- often LED pulse is different for first trial
showPhoto=1; % Show photodiode signal with LFP traces
useLED=1; % Separate trials according to LED on or off
downSampFactor=10; % Factor for down-sampling data to speed up analysis
nTrialsToUseForSpec=10; % Number of trials to use for spectrogram average
waveletDir='E:\Results\Gabor_Morlet_Wavelets\'; % For spectrogram analysis, use these saved wavelets
skipFirstTrial=1; % Sometimes LED output is wrong for the first trial - need to fix this at data acquisition at some point

% Create output directory for figures
parentSaveDir=strcat('E:\Results\NRT In Vivo Experiments\',expt.name);
if ~exist(parentSaveDir,'dir')
    mkdir(parentSaveDir);
end

% Initialize parameters for analyses
params=getAnalysisParams(expt,stimBlocks);

% For each block of constant stim. parameters, perform unit and LFP
% analyses
% Check with user to be sure all found DAQ files should be included, or
% get which should be excluded
for i=1:length(params)
    %try
        if isempty(excludeTheseDaqs)
            disp(strcat('Stim. Type:', params(i).stimType));
            disp(strcat('DAQ Files:', num2str(params(i).daqFiles)));
            exclude=input('Exclude the following DAQ files from this group. Enter answer in the format of an integer array (i.e., [1 2 3]).');
        else
            exclude=excludeTheseDaqs{i};
        end
        useTheseDaqFiles=[];
        for j=params(i).daqFiles
            if any(j==exclude)
            else
                useTheseDaqFiles=[useTheseDaqFiles j];
            end
        end
        if isempty(useTheseDaqFiles)
            continue;
        end
        params(i).daqFiles=useTheseDaqFiles;
        
        % KR hack
%         if i==1
%             params(i).daqFiles=[145 154 159];
%         end
        %disp(params(i).daqFiles);
%         if i==18
%             params(i).daqFiles=[141 148 151 161];
%         %elseif i==1
%             params(i).daqFiles=[149 152 158];
%         elseif i==20
%             params(i).daqFiles=[144 153 159];
%         elseif i==21
%             params(i).daqFiles=[145 147 150 157];
%         elseif i==22
%             params(i).daqFiles=[146 155 160];
%         elseif i==23
%             params(i).daqFiles=[154 156];
%         end
            
        
        % Analyze this stimulus parameter block
        analyzeResponseToStim(expt,params(i),strcat(parentSaveDir,'\Stim',num2str(i)),discardLastSweep,showPhoto,useLED,downSampFactor,nTrialsToUseForSpec,skipFirstTrial);
        close all;
    %catch me
    %    disp(['Problem with analyzing stim block ' num2str(i)]);
    %    disp(me);
    %end
end
    