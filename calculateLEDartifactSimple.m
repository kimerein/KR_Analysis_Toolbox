function artifact=calculateLEDartifactSimple(LFPbySweep,LFP_Fs,LED_onset,LED_offset,includeTheseTrials,plotArtifact,saveArtifact,saveName)
% Obsolete version of calculateLEDartifact.m

artifact=zeros(1,size(LFPbySweep,2));
onsetInd=floor(LED_onset*LFP_Fs)+1;
offsetInd=floor(LED_offset*LFP_Fs);
%artifact(onsetInd:offsetInd)=mean(LFPbySweep(includeTheseTrials,onsetInd:offsetInd),1);
artifact=mean(LFPbySweep(includeTheseTrials,onsetInd:offsetInd),1);

% Make initial value of artifact 0
artifact=artifact-artifact(1);

% DCtimeperiod=[1.005 1.05];
% 
% %artifact(floor(DCtimeperiod(1)*LFP_Fs)+1-onsetInd+1:floor(DCtimeperiod(2)*LFP_Fs)-onsetInd+1)=mean(artifact(floor(DCtimeperiod(1)*LFP_Fs)+1-onsetInd+1:floor(DCtimeperiod(2)*LFP_Fs)-onsetInd+1),2);
% DCvalueattime=1.005;
% artifact(floor(DCtimeperiod(1)*LFP_Fs)+1-onsetInd+1:floor(DCtimeperiod(2)*LFP_Fs)-onsetInd+1)=artifact(floor(DCvalueattime*LFP_Fs)+1-onsetInd+1)*ones(1,size(artifact(floor(DCtimeperiod(1)*LFP_Fs)+1-onsetInd+1:floor(DCtimeperiod(2)*LFP_Fs)-onsetInd+1)));
% 
% DCendInd=floor(DCtimeperiod(2)*LFP_Fs);
% x=floor((DCtimeperiod(2)+0.05)*LFP_Fs)-DCendInd;
% startY=artifact(DCendInd-onsetInd+1);
% m=startY/x;
% currVal=startY;
% for i=DCendInd+1-onsetInd+1:DCendInd+x-onsetInd+1
%     currVal=currVal-m;
%     artifact(i)=currVal;
% end
% artifact(DCendInd+x-onsetInd+2:end)=zeros(1,size(artifact(DCendInd+x-onsetInd+2:end)));

if plotArtifact==1
    figure;
    plot(artifact);
end
if saveArtifact==1
    save(saveName, 'artifact');
end