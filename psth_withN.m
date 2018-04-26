function [varargout] = psth_withN(spikes,binsize,bsmooth,duration)
% function [varargout] = psth(spikes,binsize,hAxes,bsmooth,duration)
%
% INPUTS
%   spiketimes:
%   binsize:
%   hAxes:
%   bsmooth:
%   duration: trial duration
%
% OUTPUTS
%   varargout{1} = hPsth;
%   varargout{2} = hAxes;
%   varargout{3} = n;
%   varargout{4} = centers;
%   varargout{5} = edges;
%   varargout{6} = xpoints;
%   varargout{7} = ypoints;

% Created:  3/14/10 - SRO
% Modified: 5/14/10 - SRO
%           6/8/10 - SRO
%           11/3/11 - KR passes in duration

if nargin < 2
    binsize = 50; 
end

% Use current axes if hAxes not supplied
% if nargin < 3
%     hAxes = gca;   
% end
% 
if nargin < 3
    bsmooth = 1;
end

% Set duration and number of trials
% if isfield(spikes,'sweeps')
%     a=unique(spikes.trials);
%     if length(spikes.sweeps.trials)==4*length(a)
%         numtrials=length(unique(spikes.trials));
%     else
%         numtrials = length(spikes.sweeps.trials);
%     end
% else
%     numtrials=length(unique(spikes.trials));
% end
numtrials=length(unique(spikes.trials));

% Set spiketimes
spiketimes = spikes.spiketimes;

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

% nsForStdev=zeros(length(unique(spikes.trials)),size(n,2));
% allTrials=unique(spikes.trials);
% for i=1:length(allTrials)
%     cspikes=filtspikes(spikes,0,'trials',allTrials(i));
%     nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
% end
% nsForStdev=nsForStdev/binsize;

if all(isnan(n))
    n = 0;
end

% Compute center of bins
centers = edges + diff(edges(1:2))/2;

% Last point of n contains values falling on edge(end) -- usually zero
if bsmooth
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

% Outputs
varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = numtrials;
% varargout{6} = 1*std(nsForStdev(:,1:end-1),1);