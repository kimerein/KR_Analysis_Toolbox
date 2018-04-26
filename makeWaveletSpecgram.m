function [p,LFPspecgram,allTrialSpecgrams,otherDetails]=makeWaveletSpecgram(data,trials,Fs)

global specgramOutput
% Prepare data format from KR's analysis code for use with Gabor-Morlet
% wavelet spectrogram analysis

sig=cell(length(trials),1);
for i=1:length(trials)
    sig{i}=data(trials(i),:)';
%     sig{i}=data(i,:)';
end
makeNewWavelets=1;
[p,specgramStruct,LFPspecgram,allTrialSpecgrams,otherDetails]=trySpecgram(sig,Fs,makeNewWavelets);
% save('Z:\Analysis Computer\120904_Mawake10_specgramResults\specgramStruct.mat','specgramStruct');
specgramOutput=specgramStruct;