function calcMeanAndStdevDuringWindow(spikes,window,Fs,binsize)

% Spikes must contain fields:
%   spiketimes
%   trials
% window(1) is start of window for which to calculate mean and stdev of FR
% across trials; window(2) is end of window
% Fs is sampling rate of the spikes data

[n,centers,edges,xpoints,ypoints,stds]=psth_wStdev_valuesOnly(spikes,2);
m=mean(ypoints(xpoints>window(1) & xpoints<window(2)));
s=stdev(