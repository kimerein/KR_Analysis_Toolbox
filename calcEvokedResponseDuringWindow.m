function [varargout]=calcEvokedResponseDuringWindow(spikes,window,baseline_window)
%function [varargout] = psth_wStdev_valuesOnly(spikes,binsize,bsmooth)

% Spikes must contain fields:
%   spiketimes
%   trials
% window(1) is start of window for which to calculate mean and stdev of FR
% across trials; window(2) is end of window
% Fs is sampling rate of the spikes data

binsize=1; % in ms

% Set duration and number of trials
% duration = 6;     % Maybe add checking for equal sweep durations?
if isfield(spikes,'sweeps')
    numtrials = length(spikes.sweeps.trials);
else
    numtrials=length(unique(spikes.trials));
end

% Set spiketimes
spiketimes = spikes.spiketimes;

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
%edges = 0:binsize:duration;
edges=window(1):binsize:window(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

edges_base=baseline_window(1):binsize:baseline_window(2);
n_base = histc(spiketimes,edges_base);
n_base = n_base/numtrials/binsize;
m_base = mean(n_base);

nsForStdev=zeros(length(unique(spikes.trials)),size(n,2));
nsForStdev_base=zeros(length(unique(spikes.trials)),size(n_base,2));
allTrials=unique(spikes.trials);
for i=1:length(allTrials)
    cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
    nsForStdev_base(i,:)=histc(cspikes.spiketimes,edges_base);
end
nsForStdev=nsForStdev/binsize;
% nsForStdev=[nsForStdev; zeros(numtrials-length(allTrials),size(nsForStdev,2))]/binsize;
nsForStdev_base=nsForStdev_base/binsize;
% nsForStdev_base=[nsForStdev_base; zeros(numtrials-length(allTrials),size(nsForStdev_base,2))]/binsize;
s=std(mean(nsForStdev,2),1);

if all(isnan(n))
    n = 0;
end

varargout{1} = m-m_base;
varargout{2} = s;
varargout{3} = mean(nsForStdev,2)-mean(nsForStdev_base,2);