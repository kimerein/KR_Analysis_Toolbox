function fig=makeSingleSpectrogram(data,Fs,useTheseTrials,makeNewWavelets)

disp('Making spectrogram');
if isempty(useTheseTrials)
   return
end
data=data(useTheseTrials,:);
sig=cell(size(data,1),1);
for i=1:size(data,1)
    sig{i}=data(i,:)';
end
[fig,specgramStruct]=trySpecgram(sig,Fs,makeNewWavelets);
save(strcat('E:\Results\GammaLFPs\Spectrograms\specgramStruct13.mat'),'specgramStruct');
