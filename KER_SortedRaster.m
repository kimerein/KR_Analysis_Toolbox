function [varargout] = KER_SortedRaster(spikes,redTrials,nParams1,nParams2,trialLength,black)
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
spacing = 0.03;

% Make raster
hRaster = zeros(size(spiketimes));

hold off;

lineColors=zeros(100,3);
l=1;
for m=0:0.16:0.8
    for j=0:0.16:0.8
        for k=0:0.16:0.8
            lineColors(l,:)=[m j k];
            l=l+1;
        end
    end
end
lineColors=lineColors(randperm(size(lineColors,1)),:);
lineColorN=1;

mostTrialsForStim=0;
for i=1:nParams1
    for j=1:nParams2
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*nParams1+j);
        trialsForStim=length(unique(someSpikes.trials));
        if trialsForStim>mostTrialsForStim
            mostTrialsForStim=trialsForStim;
        end
    end
end

countStims=0;
countRows=0;
countSpikes=1;
for i=1:nParams1
    for j=1:nParams2
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*nParams1+j);
        %(i-1)*nParams1+j
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
                    %(countStims*nParams2)+cTrials
                end
                x = [spiketimes(k) spiketimes(k)];
                %y = [(countStims*nParams2)+cTrials-spacing ((countStims*nParams2)+cTrials-1+spacing)];
                y = [countRows*spacing countRows*spacing+spacing];
                if black==1
                    hRaster(countSpikes) = line(x,y,'Color','k','LineWidth',1);
                else
                    hRaster(countSpikes) = line(x,y,'Color',lineColors((i-1)*nParams1+j,:),'LineWidth',1);
                end
                countSpikes=countSpikes+1;
            end
        end  
        if rowsForThisStim<mostTrialsForStim
            countRows=countRows+(mostTrialsForStim-rowsForThisStim);
        end
        countStims=countStims+1;
    end
    %line([0 trialLength],[((countStims*nParams2)+cTrials+spacing) ((countStims*nParams2)+cTrials+spacing)],'Color','k');
    line([0 trialLength],[countRows*spacing+spacing countRows*spacing+spacing],'Color','k');
end

% for i = 1:length(spiketimes)
%     x = [spiketimes(i) spiketimes(i)];
%     y = [trials(i)-spacing (trials(i)-1+spacing)];
%     if any(trials(i)==redTrials)
%         c='r';
%     else
%         c='b';
%     end
% %     hRaster(i) = line(x,y,'Color',c,'Parent',hAxes,'LineWidth',1);
%     hRaster(i) = line(x,y,'Color',c,'LineWidth',1);
% end

% Set properties
%numtrials = max(trials);
%numtrials = (countStims*nParams2)+cTrials;
numtrials=countRows;
offset = numtrials*0.03;
% set(hAxes,'TickDir','out','YDir','reverse','FontSize',9, ... 
%     'YLim',[(0-offset) (numtrials+offset)]);
set(hAxes,'TickDir','out','YDir','reverse','FontSize',9, ... 
    'YLim',[0 offset]);
xlabel(hAxes,'seconds')
ylabel(hAxes,'trials')

% Outputs
varargout{1} = hAxes;
varargout{2} = hRaster;

