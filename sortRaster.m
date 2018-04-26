function [varargout] = sortRaster(spikes,params,black,colorCode,fig)
% INPUTS
%   spikes: The spikes struct. This function requires that spikes contain
%   only the following fields:
%   .spiketimes: Vector of spike times
%   .trials: Vector indicating trial in which each spike occurred
%
%   params are the relevant stimulus parameters
%   must contain:
%       params.Var1_values
%       params.Var2_values
%
%   if black is 1, all ticks will be black
%   else if black is 0, the colorCode will be used to decide the color for
%   each row
% 
%   colorCode is a structure array containing the fields:
%       colorCode.color = a three-element vector (values between 0 and 1)
%                         describing the color associated with the below
%                         stimulus conditions
%       colorCode.var1_values = an array of Stim. Var. 1 values which (in
%                               conjunction with var2_values) are represented
%                               by the color 
%       colorCode.var2_values = an array of Stim. Var. 2 values which (in
%                               conjunction with var1_values) are represented
%                               by the color 
%   (see makeColorCode)
%
% OUTPUTS
%   varargout(1) = hAxes
%   varargout(2) = hRaster

% Set spiketimes and trial number
spiketimes = spikes.spiketimes;

% Get current axes
set(0,'CurrentFigure',fig);
hAxes = gca;

% Set raster line properties
spacing = 0.03;

% Make raster
hRaster = zeros(size(spiketimes));

hold off;

% Find the max iterations for a stimulus condition
% For stimulus conditions that do not have that many iterations, 
% add blank rows to raster
mostTrialsForStim=0;
for i=1:length(params.Var1_values)
    for j=1:length(params.Var2_values)
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*params.Var1_values+j);
        trialsForStim=length(unique(someSpikes.trials));
        if trialsForStim>mostTrialsForStim
            mostTrialsForStim=trialsForStim;
        end
    end
end

% Make raster
countStims=0;
countRows=0;
countSpikes=1;
for i=1:length(params.Var1_values)
    for j=1:length(params.Var2_values)
        currColor=colorCode((i-1)*length(params.Var1_values)+j).color;
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*length(params.Var1_values)+j);
        spiketimes = someSpikes.spiketimes;
        rowsForThisStim=0;
        if ~isempty(spiketimes)
            rowsForThisStim=1;
            countRows=countRows+1;
            currTrial=min(someSpikes.trials);
            cTrials=1;
            for k=1:length(spiketimes)
                if someSpikes.trials(k)~=currTrial
                    cTrials=cTrials+1;
                    currTrial=someSpikes.trials(k);
                    countRows=countRows+1;
                    rowsForThisStim=rowsForThisStim+1;
                end
                x = [spiketimes(k) spiketimes(k)];
                y = [countRows*spacing countRows*spacing+spacing];
                if black==1
                    hRaster(countSpikes) = line(x,y,'Color','k','LineWidth',1);
                else
                    hRaster(countSpikes) = line(x,y,'Color',currColor,'LineWidth',1);
                end
                countSpikes=countSpikes+1;
            end
        end  
        if rowsForThisStim<mostTrialsForStim
            countRows=countRows+(mostTrialsForStim-rowsForThisStim);
        end
        countStims=countStims+1;
    end
    line([0 params.totalTrialLength],[countRows*spacing+spacing countRows*spacing+spacing],'Color','k');
end

% Set properties
numtrials=countRows;
offset = numtrials*0.03;
set(hAxes,'TickDir','out','YDir','reverse','FontSize',9, ... 
    'YLim',[0 offset]);
xlabel(hAxes,'seconds')
ylabel(hAxes,'trials')

% Outputs
% varargout{1} = hAxes;
% varargout{2} = hRaster;
varargout{1}=fig;

