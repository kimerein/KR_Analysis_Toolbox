

disp('Making spectrogram');
data=m';
sig=cell(size(data,1),1);
for i=1:size(data,1)
    sig{i}=data(i,:)';
end
[fig,specgramStruct]=trySpecgram(sig,Fs,makeNewWavelets);
save(strcat('E:\Results\GammaLFPs\Spectrograms\specgramStruct4.mat'),'specgramStruct');
