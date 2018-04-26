function [varargout] = old_psth(spikes,binsize,duration)
% function [varargout] = psth(spiketimes,binsize,duration)
% INPUTS
%   spiketimes:
%   binsize:
%   duration:
%
% OUTPUTS
%   varargout{1} = hAxes:
%   varargout{2} = hPsth:
%   varargout{3} = n:
%   varargout{4} = edges:
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
    spikes = spikes.spiketimes;% BA switched to spikes
end

% Get current axes
hAxes = gca;

% Convert binsize from s to ms
binsize = binsize/1000;

% Make PSTH
edges = 0:binsize:duration;
[n,bin] = histc(spikes,edges);
n = n/length(trials)/binsize;
% n = n/12/binsize;
hPsth = line(edges,n,'Parent',hAxes,'LineWidth',1.5);  %temp

% Outputs
varargout{1} = hAxes;
varargout{2} = hPsth;
varargout{3} = n;
varargout{4} = edges;