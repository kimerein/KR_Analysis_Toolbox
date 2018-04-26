function [varargout] = K_psth(spikes,binsize,duration,lineColor)
% INPUTS
%   spiketimes: a field of spikes structure with times of spikes
%   binsize: size of time bin in seconds for which to calculate firing rate
%   duration: trial duration
%
% OUTPUTS
%   varargout{1} = hAxes: handle to the axes of the psth
%   varargout{2} = hPsth: handle to the psth
%   varargout{3} = n: the firing rate in each bin
%   varargout{4} = edges: the times the separate bins
%
% 3/14/10 - SRO

if nargin < 3
    duration = max(spikes.spiketimes);
end
if nargin < 2
    binsize = 50; 
end

trials = spikes.sweeps.trials;

% Set spiketimes
if isstruct(spikes)
    spikes = spikes.spiketimes;
end

% Get current axes
hAxes = gca;

% Make PSTH
edges = 0:binsize:duration;
[n,bin] = histc(spikes,edges);
n = n/length(trials)/binsize;
hPsth = line(edges,n,'Parent',hAxes,'LineWidth',1.5,'Color',lineColor);  %temp

% Outputs
varargout{1} = hAxes;
varargout{2} = hPsth;
varargout{3} = n;
varargout{4} = edges;