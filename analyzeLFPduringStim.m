function [saveToDir,analysisInfo]=analyzeLFPduringStim(expt,daqs1,daqs2,compare,params1,params1Name,params2,params2Name,currSaveDir,extraN,block1Trials,block2Trials)
% Analyze the LFP responses to stimuli
% compare is 0 or 1
% 0 if there should NOT be a comparison of LFP block1Trials and
% block2Trials, 1 if the analysis should compare these
%
% if daqs1 is not empty, read in these daqNames and concatenate the data
% else use the daq file numbers in params1
% same for daqs2
%
% block1Trials contains indices into the params1 data
%
% block2Trials contains indices into the params2 data
%
% params1 are the stimulus parameters for the first type of stimulus;
% same for params2
%
% COMPARED LFP BLOCKS MUST HAVE SAME SAMPLING RATE FOR THIS CODE TO WORK
% add a fix later for this


global physChannel photoChannel LFPdef_lowerCutoff LFPdef_upperCutoff BPLowerCutoff BPUpperCutoff;

% Get DAQ files to include in the analysis
if isempty(daqs1)
    daqFileNames1=expt.files.names(params1.daqFiles);
else
    daqFileNames1=daqs1;
end
if compare==1
    if isempty(daqs2)
        daqFileNames2=expt.files.names(params2.daqFiles);
    else
        daqFileNames2=daqs2;
    end
end

% Make output directory for LFP
if compare==1
    saveToDir=[currSaveDir '\' params1Name '_vs_' params2Name num2str(extraN)];
else
    saveToDir=strcat(currSaveDir,'\',params1Name,num2str(extraN));
end
if ~exist(saveToDir,'dir')
    mkdir(saveToDir);
end

% See if photodiode-aligned data already exists
% If so, load
% Else get data and align to photodiode
photoData1=what(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name));
if compare==1
    photoData2=what(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name));
end
% Directory should contain photodiode-aligned LFPbySweep,
% bandPassedLFPbySweep, photodiodeBySweep
if isempty(photoData1)
    alignData1=1;
elseif isempty(photoData1.mat)
    alignData1=1;
elseif length(photoData1.mat)<4
    alignData1=1;
else
    alignData1=0;
end
if compare==1
    if isempty(photoData2)
        alignData2=1;
    elseif isempty(photoData2.mat)
        alignData2=1;
    elseif length(photoData2.mat)<4
        alignData2=1;
    else
        alignData2=0;
    end
end
if alignData1==1   
    % Get LFP 
    [LFPbySweep1,LFP_Fs1,bandPassedLFPbySweep1,photodiodeBySweep1]=getLFPbySweepSeveralFiles(daqFileNames1,params1.samplingRate,LFPdef_lowerCutoff,LFPdef_upperCutoff,BPLowerCutoff,BPUpperCutoff,physChannel,photoChannel);

    % Align to photodiode
    [LFPbySweep1,LFP_Fs1,bandPassedLFPbySweep1,photodiodeBySweep1,params1]=alignLFPToPhotodiode(LFPbySweep1,LFP_Fs1,bandPassedLFPbySweep1,photodiodeBySweep1,params1);

    % Save data aligned to photodiode
    photoAlignedDir=strcat(currSaveDir,'\PhotoAligned_Data\',params1Name);
    if ~exist(photoAlignedDir,'dir')
        mkdir(photoAlignedDir);
    end
    save(strcat(photoAlignedDir,'\LFPbySweep.mat'),'LFPbySweep1');
    save(strcat(photoAlignedDir,'\bandPassedLFPbySweep.mat'),'bandPassedLFPbySweep1');
    save(strcat(photoAlignedDir,'\photodiodeBySweep.mat'),'photodiodeBySweep1');
    save(strcat(photoAlignedDir,'\LFP_Fs.mat'),'LFP_Fs1');
    p=params1.ONstart;
    save(strcat(photoAlignedDir,'\ONstart.mat'),'p');
else
    % Just read in mat files
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name,'\LFPbySweep.mat'));
    LFPbySweep1=a.LFPbySweep1;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name,'\bandPassedLFPbySweep.mat'));
    bandPassedLFPbySweep1=a.bandPassedLFPbySweep1;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name,'\photodiodeBySweep.mat'));
    photodiodeBySweep1=a.photodiodeBySweep1;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name,'\LFP_Fs.mat'));
    LFP_Fs1=a.LFP_Fs1;
    a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params1Name,'\ONstart.mat'));
    params1.ONstart=a.p;
end  
if compare==1
    if alignData2==1   
        % Get LFP
        [LFPbySweep2,LFP_Fs2,bandPassedLFPbySweep2,photodiodeBySweep2]=getLFPbySweep(daqFileNames2,params2.samplingRate,LFPdef_lowerCutoff,LFPdef_upperCutoff,BPLowerCutoff,BPUpperCutoff,physChannel,photoChannel);
        
        % Align to photodiode
        %[LFPbySweep2,LFP_Fs2,bandPassedLFPbySweep2,photodiodeBySweep2,params2]=alignLFPToPhotodiode(LFPbySweep2,LFP_Fs2,bandPassedLFPbySweep2,photodiodeBySweep2,params2);
        
        % Save data aligned to photodiode
        photoAlignedDir=strcat(currSaveDir,'\PhotoAligned_Data\',params2Name);
        if ~exist(photoAlignedDir,'dir')
            mkdir(photoAlignedDir);
        end
        save(strcat(photoAlignedDir,'\LFPbySweep.mat'),'LFPbySweep2');
        save(strcat(photoAlignedDir,'\bandPassedLFPbySweep.mat'),'bandPassedLFPbySweep2');
        save(strcat(photoAlignedDir,'\photodiodeBySweep.mat'),'photodiodeBySweep2');
        save(strcat(photoAlignedDir,'\LFP_Fs.mat'),'LFP_Fs2');
        p=params1.ONstart;
        save(strcat(photoAlignedDir,'\ONstart.mat'),'p');
    else
        % Just read in mat files
        a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name,'LFPbySweep.mat'));
        LFPbySweep2=a.LFPbySweep2;
        a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name,'bandPassedLFPbySweep.mat'));
        bandPassedLFPbySweep2=a.bandPassedLFPbySweep2;
        a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name,'photodiodeBySweep.mat'));
        photodiodeBySweep2=a.photodiodeBySweep2;
        a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name,'LFP_Fs.mat'));
        LFP_Fs2=a.LFP_Fs2;
        a=load(strcat(currSaveDir,'\PhotoAligned_Data\',params2Name,'ONstart.mat'));
        params2.ONstart=a.p;
    end
end

% Make sure sampling rates match
% Write code later
LFP_Fs=LFP_Fs1;

% Put block1 and block2 trials together, cutting the lengths of each sweep down
% to the minimum size for block1 and block2 sweeps
if ~isempty(block1Trials)
    if block1Trials==-100
        block1Trials=1:size(LFPbySweep1,1);
    end
end
if ~isempty(block2Trials)
    if block2Trials==-100
        block2Trials=1:size(LFPbySweep2,1);
    end
end
if compare==1
    lengthBlock1=size(LFPbySweep1,2);
    lengthBlock2=size(LFPbySweep2,2);
    blockLength=min(lengthBlock1,lengthBlock2);    
    LFPbySweep=[LFPbySweep1(:,1:blockLength); LFPbySweep2(:,1:blockLength)];
    bandPassedLFPbySweep=[bandPassedLFPbySweep1(:,1:blockLength); bandPassedLFPbySweep2(:,1:blockLength)];
else
    LFPbySweep=LFPbySweep1;
    bandPassedLFPbySweep=bandPassedLFPbySweep1;
end
    
 % LFP Analyses
[LFPfigs,usedRed,usedBlue]=analyzeLFPSweepsForRinging(compare,LFPbySweep,bandPassedLFPbySweep,LFP_Fs,block1Trials,block2Trials+length(block1Trials),params1); 
analysisInfo.usedRed=usedRed;
analysisInfo.usedBlue=usedBlue-length(block1Trials);
% Save figures from this analysis and the analysis details used to
% produce them
saveas(LFPfigs(1),strcat(currSaveDir,'\','stimulusTriggeredLFP.bmp')); 
saveas(LFPfigs(1),strcat(currSaveDir,'\','stimulusTriggeredLFP.fig'));  
save(strcat(currSaveDir,'\','AnalysisDetails.mat'),'analysisInfo');

% sig=cell(1,1);
% sig{1}=mean(LFPbySweep,1)';
% trySpecgram(sig,LFP_Fs);
% 'hi'