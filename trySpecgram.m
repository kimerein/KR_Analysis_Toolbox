function [p,specgramStruct,LFPspecgram,allTrialSpecgrams,otherDetails]=trySpecgram(sig,lowPass_samplingRate,makeNewWavelets)
global waveletDir

%try
    % X_div gives the time between tick marks for the spectrogram
    X_div=0.5; % in seconds
    % Y_div gives the frequency bands between tick marks for the
    % spectrogram
    Y_div=10; % in Hz
    if exist([waveletDir 'gbwavelets_freq.mat'],'file') && exist([waveletDir 'gbwavelets_gabor.mat'],'file') && makeNewWavelets~=1
        disp('Using existing wavelets.');
        [gb_params, freq, gabor]=gabor_morlet_config(lowPass_samplingRate,[waveletDir 'gbwavelets_freq.mat'],[waveletDir 'gbwavelets_gabor.mat']);
        %[gb_params, freq, gabor]=gabor_morlet_config(1000,[waveletDir 'gbwavelets_freq.mat'],[waveletDir 'gbwavelets_gabor.mat']);
    else
        disp('Making new wavelets.');
        [gb_params, freq, gabor]=gabor_morlet_config(lowPass_samplingRate,[],[]);
        disp('Finished making new wavelets.');
        %[gb_params, freq, gabor]=gabor_morlet_config(1000,[],[]);
    end
    if length(sig)==1
        if any(size(sig{1})==0)
            p=[];
            return
        end
        %p=gabor_morlet_plot(sig,gb_params,freq,gabor,1,sig{1},'LFP',X_div,Y_div,1,[1 10 20 30 40 50 60 80 100],lowPass_samplingRate);
        allTrialSpecgrams=cell(1,1);
        [p,specgramStruct,LFPspecgram]=gabor_morlet_plot(sig,gb_params,freq,gabor,0,sig{1},'LFP',X_div,Y_div,1,[1 10 20 30 40 50 60 80 100],lowPass_samplingRate);
        allTrialSpecgrams{1}=LFPspecgram;
%         disp('hi');
    else
        [p,specgramStruct,LFPspecgram,allTrialSpecgrams]=gabor_morlet_plot(sig,gb_params,freq,gabor,0,[],'LFP',X_div,Y_div,1,[1 10 20 30 40 50 60 80 100],lowPass_samplingRate);
    end
%catch me
%     disp('Error making spectrogram');
%     p=[];
%end

otherDetails.gb_params=gb_params;
otherDetails.freq=freq;
otherDetails.gabor=gabor;