function OFFperiod_data=KER_analyzeLFPringing
daqFileNames=cell(1,2);
daqFileNames{1}='KR_2010-05-14_M3D343_1.daq';
daqFileNames{2}='KR_2010-05-14_M3D343_2.daq';
cellName='Mouse 3 LFP Stim 1';
ONlength=0.5;
totalTrialLength=1.5;
takeNSweeps=35;
saveToDir=strcat('C:\Documents and Settings\Admin\My Documents\Recurrent Activity - Data\Data from Expt. on 05-14-10\',cellName);
if ~exist(saveToDir,'dir')
    mkdir(saveToDir);
end

if ~exist('LFPbySweep','var')
    LFPbySweep=[];
end

[LFPbySweep,Fs,bandPassedLFPbySweep]=KER_analyzeLFP_forRinging(daqFileNames,LFPbySweep,ONlength,totalTrialLength,takeNSweeps);

saveToDir=strcat(saveToDir,'\');
saveas(figure(2),strcat(saveToDir,'nEquals35_stimTriggeredLFP_withPhotodiode.bmp'));
saveas(figure(3),strcat(saveToDir,'nEquals35_LFP_stimTriggeredLFP.bmp'));
saveas(figure(4),strcat(saveToDir,'nEquals35_LFP_stimTriggeredLFP_withStdDev.bmp'));
saveas(figure(6),strcat(saveToDir,'nEquals35_LFP_stimTriggeredLFP_filtered30to80Hz.bmp'));
saveas(figure(8),strcat(saveToDir,'nEquals35_LFP_stimTriggeredLFP_filtered30to80Hz_withStdDev.bmp'));
%saveas(figure(7),strcat(saveToDir,'stimTriggeredLFP_filtered30to80Hz_ampDiffWithGreaterThanEqualTo',num2str(redThresh),'OFFSpikes.bmp'));


figure;
plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep,1),'Color','k');
axis([ONlength ONlength+1 -0.02 0.02]);

saveas(figure(10),strcat(saveToDir,'nEquals35_stimTriggeredLFP_OFFgamma.bmp'));

disp('Number of Trials');
size(LFPbySweep,1)

end