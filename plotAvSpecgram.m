function [p,specgramStruct,LFPspecgram]=plotAvSpecgram(trialSpecgrams,Nsamples,sampling_rate,freq)

X_div=0.5; % in seconds
Y_div=10; % in Hz
yes_logY=1;
Y_ticks=[1 10 20 30 40 50 60 80 100];
for i=1:length(trialSpecgrams.specgrams)
    if i==1
        tf=trialSpecgrams.specgrams{i};
    else
        tf=tf+trialSpecgrams.specgrams{i};
    end
end
tf=tf/length(trialSpecgrams.specgrams);

mtf=repmat(mean(tf,2),1,size(tf,2));
xtf=(tf-mtf)./(tf+mtf);

figure(); 
plot_time_freq(tf,0,Nsamples/sampling_rate,X_div,freq(1),freq(end),Y_div,yes_logY,Y_ticks);
ti=sprintf('Time-Frequency Analysis of Data');
title(ti);

p=figure(); 
plot_time_freq(xtf,0,Nsamples/sampling_rate,X_div,freq(1),freq(end),Y_div,yes_logY,Y_ticks);

specgramStruct.plot_time_freq_arg1=xtf;
specgramStruct.plot_time_freq_arg2=0;
specgramStruct.plot_time_freq_arg3=Nsamples/sampling_rate;
specgramStruct.plot_time_freq_arg4=X_div;
specgramStruct.plot_time_freq_arg5=freq(1);
specgramStruct.plot_time_freq_arg6=freq(end);
specgramStruct.plot_time_freq_arg7=Y_div;
specgramStruct.plot_time_freq_arg8=yes_logY;
specgramStruct.plot_time_freq_arg9=Y_ticks;
LFPspecgram=tf;