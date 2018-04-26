function [LFPbySweep,Fs,bandPassedLFPbySweep,photoData]=getLFPbySweepSeveralFiles(daqFileNames,Fs,def_LowerCutoff,def_UpperCutoff,lower_cutoff,upper_cutoff,physChannel,photoChannel,ledChannel)
% Reads in several daq files and concatenates data
% GLOBAL
% dataDir is a string specifying the directory containing the data to be
% read in
% INPUT PARAMETERS
% daqFileNames              a cell array of strings giving the names of the
%                           .daq files to be read in
% Fs                        sampling rate of the daq files
% def_LowerCutoff           all data read in from the daq files will be
%                           band-pass filtered between def_LowerCutoff and
%                           def_UpperCutoff
% def_UpperCutoff           "
% lower_cutoff              the returned array bandPassedLFPbySweep will be
%                           the data from the daq files sorted into trials and 
%                           band-pass filtered between lower_cutoff and upper_cutoff
% physChannel               an index specifying the physiology channel in
%                           the daq files (was used when the daq file was created)
% photoChannel              an index specifying the photodiode channel in
%                           the daq files 
% ledChannel                an index specifying the LED channel in the daq
%                           files
% OUTPUT
% LFPbySweep                a matrix of the data from the daq files (only
%                           including in the physChannel and photoChannel
%                           channels in the daq files), sorted into sweeps
%                           each row is a trial/sweep, each column is a
%                           sample
% Fs                        sampling rate of LFPbySweep and bandPassedLFPbySweep and photoData
%                           after daq file data has been down-sampled (w/
%                           down sampling factor of 10)
% bandPassedLFPbySweep      LFPbySweep data but filtered between
%                           lower_cutoff and upper_cutoff
% photoData                 a matrix of photodiode sorted into sweeps

global dataDir

LFPbySweep=[];
photoData=[];
for i=1:length(daqFileNames)
    daqFileName=daqFileNames{i};
    [thisLFPbySweep,Fs2,thisPhotoData]=readStimTriggeredLFP(daqFileName,Fs,0,def_LowerCutoff,def_UpperCutoff,10,physChannel,photoChannel,ledChannel);
    if i>1
        % Truncate these trials to match the first loaded file
        if size(thisLFPbySweep,2)>size(LFPbySweep,2)
            thisLFPbySweep=thisLFPbySweep(:,1:size(LFPbySweep,2));
            thisPhotoData=thisPhotoData(:,1:size(thisPhotoData,2));
        end
    end
    LFPbySweep=[LFPbySweep; thisLFPbySweep];
    photoData=[photoData; thisPhotoData];
end
Fs=Fs2;

bandPassedLFPbySweep=bandPassLFP(LFPbySweep,Fs,lower_cutoff,upper_cutoff,1);
end