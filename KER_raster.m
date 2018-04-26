function [varargout] = KER_raster(spikes,redTrials,blueTrials)
% INPUTS
%   spikes: The spikes struct. This function requires that spikes contain
%   only the following fields:
%   .spiketimes: Vector of spike times
%   .trials: Vector indicating trial in which each spike occurred
%   redTrials are the set of trials to be colored red
%   all other trials will be colored blue
%
% OUTPUTS
%   varargout(1) = hAxes
%   varargout(2) = hRaster

% Set spiketimes and trial number
spiketimes = spikes.spiketimes;
trials = spikes.trials;

% Get current axes
figure;
hAxes = gca;

% Set raster line properties
spacing = 0.02;

% Make raster
hRaster = zeros(size(spiketimes));

hold off;

for i = 1:length(spiketimes)
    x = [spiketimes(i) spiketimes(i)];
    if spiketimes(i)>3.5
        disp(trials(i));
    end
    y = [trials(i)-spacing (trials(i)-1+spacing)];
    if any(trials(i)==redTrials)
        c='r';
    else
        c='b';
    end
%     hRaster(i) = line(x,y,'Color',c,'Parent',hAxes,'LineWidth',1);
    hRaster(i) = line(x,y,'Color',c,'LineWidth',1);
end

% Set properties
% numtrials = max(trials);
numtrials = length(spikes.sweeps.trials);
offset = numtrials*0.03;
set(hAxes,'TickDir','out','YDir','reverse','FontSize',9, ... 
    'YLim',[(0-offset) (numtrials+offset)])
xlabel(hAxes,'seconds')
ylabel(hAxes,'trials')

% Outputs
varargout{1} = hAxes;
varargout{2} = hRaster;

