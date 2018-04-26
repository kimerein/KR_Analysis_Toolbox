function [cycleAvs_con,cycleAvs_led,stimOverBase_con,stimOverBase_led]=plotCycleAverage(expt,spikes,useAssign,nBins,closeFig)
% Spikes should be all spikes so that you can count all trials

if ~isempty(useAssign)
    spikes=filtspikes(spikes,0,'assigns',useAssign);
end

stimFreqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
amberFreqs=[1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
% stimFreqs=[1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05 16.05 18.05 20.05 30.05 40.05 50.05 60.05];
% amberFreqs=[1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05 16.05 18.05 20.05 30.05 40.05 50.05 60.05];

stimWindow=[1 4];
baseWindow=[0 1];

getStimCycles=cell(1,length(stimFreqs));
getNullCycles=cell(1,length(stimFreqs));
for i=1:length(stimFreqs)
    % stimFreqs in cycles per second
    lengthCycle=1/stimFreqs(i); % in s
    getStimCycles{i}=stimWindow(1):lengthCycle:stimWindow(2);
    getNullCycles{i}=baseWindow(1):lengthCycle:baseWindow(2);
end

% Figure layout
h.fig = landscapeFigSetup;
set(h.fig,'Visible','off','Position',[792 399 1056 724]);

useNcycles=-1;
% Get spiketimes and trials for each stimulus frequency
cycleAvs_con=cell(1,length(stimFreqs));
cycleAvs_led=cell(1,length(amberFreqs));
stimOverBase_con=zeros(1,length(stimFreqs));
stimOverBase_led=zeros(1,length(stimFreqs));
for i=1:length(stimFreqs)
    h.r.ax(i)=axes;
    h.psth.ax(i)=axes;
    defaultAxes(h.r.ax(i));
    temp=getNullCycles{i};
    nwait=temp(2)-temp(1);
    temp=getStimCycles{i};
    nstim=temp(2)-temp(1);
    dur=nwait+nstim;
    % Get spiketimes and trials
    % Control trials
    % psth function gets numtrials count from length(spikes.sweeps.trials)
    cspikes=makeTempField(spikes,'led',stimFreqs(i));
    cspikes=filtspikes(cspikes,0,'temp',1);
    for j=2:length(getNullCycles{i})-1
        if useNcycles~=-1 && j>=useNcycles+1
            n=getNullCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=-1;
        else
            n=getNullCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))-n(j);
        end
    end
    a=unique(cspikes.spiketimes);
    cspikes=filtspikes(cspikes,0,'spiketimes',a(a~=-1));
    for j=2:length(getStimCycles{i})-1
        if useNcycles~=-1 && j>=useNcycles+1
            n=getStimCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=-1;
        else
            n=getStimCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))-n(j)+nwait;
        end
    end
    a=unique(cspikes.spiketimes);
    cspikes=filtspikes(cspikes,0,'spiketimes',a(a~=-1));
    h.r.l(i,1)=raster(cspikes,h.r.ax(i),1,0,dur);
    set(h.r.l(i,1),'Color','k');
    if (nstim/nBins)*1000<1
        psthBin=1;
    else
        psthBin=(nstim/nBins)*1000;
    end
    [h.psth.l(i,1),~,~,~,~,cax,cay]=psthForCycle(cspikes,psthBin,h.psth.ax(i),0,dur);
    cycleAvs_con{i}=[cax; cay];
    halfN=ceil(length(cay)/2);
    stimOverBase_con(i)=mean(cay(halfN+1:end))-mean(cay(1:halfN));
    set(h.psth.l(i,1),'Color','k');
    % Amber trials
    cspikes=makeTempField(spikes,'led',amberFreqs(i));
    cspikes=filtspikes(cspikes,0,'temp',1);
    for j=2:length(getNullCycles{i})-1
        if useNcycles~=-1 && j>=useNcycles+1
            n=getNullCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=-1;
        else
            n=getNullCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))-n(j);
        end
    end
    if useNcycles~=-1
        a=unique(cspikes.spiketimes);
        cspikes=filtspikes(cspikes,0,'spiketimes',a(a~=-1));
    end
    for j=2:length(getStimCycles{i})-1
        if useNcycles~=-1 && j>=useNcycles+1
            n=getStimCycles{i};
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=-1;
        else
            n=getStimCycles{i}; 
            cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))=cspikes.spiketimes(cspikes.spiketimes>n(j) & cspikes.spiketimes<n(j+1))-n(j)+nwait;
        end
    end
    if useNcycles~=-1
        a=unique(cspikes.spiketimes);
        cspikes=filtspikes(cspikes,0,'spiketimes',a(a~=-1));
    end
    h.r.l(i,2)=raster(cspikes,h.r.ax(i),1,0,dur);
    set(h.r.l(i,2),'Color','r');
    if (nstim/nBins)*1000<1
        psthBin=1;
    else
        psthBin=(nstim/nBins)*1000;
    end
    [h.psth.l(i,2),~,~,~,~,cax,cay]=psthForCycle(cspikes,psthBin,h.psth.ax(i),0,dur);
    cycleAvs_led{i}=[cax; cay];
    halfN=ceil(length(cay)/2);
    stimOverBase_led(i)=mean(cay(halfN+1:end))-mean(cay(1:halfN));
    set(h.psth.l(i,2),'Color','r');
    xlim(h.psth.ax(i),[0 dur]);
end
h.r.ax=h.r.ax';
h.psth.ax=h.psth.ax';
set(h.r.ax,'Box','on');
% Set axes properties
hTemp = reshape(h.r.ax,numel(h.r.ax),1);
ymax = setSameYmax(hTemp);
removeAxesLabels(hTemp);
defaultAxes(hTemp);
gray = [0.85 0.85 0.85];
set(hTemp,'YColor',gray,'XColor',gray,'XTick',[],'YTick',[]);
setRasterPSTHpos(h);
hTemp=reshape(h.psth.ax,numel(h.psth.ax),1);
% set(h.fig,'Visible','on');
ymax=setSameYmax(hTemp,15);
removeInd=1:length(hTemp);
keepInd=ceil(length(hTemp)/2) + 1;
b.xl = 1;
b.yl = 1;
b.xtl = 0;
b.ytl = 1;
removeAxesLabels(hTemp(setdiff(removeInd,keepInd)),b);
defaultAxes(hTemp,0.25,0.25);

% Compute average waveform
[h.avgwv.l h.avgwv.ax maxch] = plotAvgWaveform(spikes,0);
defaultAxes(h.avgwv.ax,0.22,0.24);
setSameYmax(h.avgwv.ax,2,1);
xlabel('ms');
ylabel('mV');
% Make autocorrelation plot
h.autocorr.ax = axes;
plotAutoCorr(spikes,h.autocorr.ax,50,1);
defaultAxes(h.autocorr.ax,0.22,0.2);

% Define locations in axes matrix
h.mat(1).params.matpos = [0 0.68 0.49 0.35];                % [left top width height]
h.mat(1).params.figmargin = [0.00 0 0 0.05];                % [left right top bottom]
h.mat(1).params.matmargin = [0 0 0 0];                      % [left right top bottom]
h.mat(1).params.cellmargin = [0.05 0.035 0.05 0.05];        % [left right top bottom]
h.mat(1).ncol = 2;
h.mat(1).nrow = 1;
h.mat(1).h(1) = h.autocorr.ax;
h.mat(1).h(2) = h.avgwv.ax;

for i=1:length(h.mat)
    ind = 1:length(h.mat(i).h);
    setaxesOnaxesmatrix(h.mat(i).h,h.mat(i).nrow,h.mat(i).ncol,ind, ...
        h.mat(i).params,h.fig);
end

% Add expt name and unit info.
exptInfo = ['T' num2str(spikes.info.trodeInd) '_' num2str(useAssign)];
h.textbox = annotation('textbox',[0 0 0.3 0.022],'String',exptInfo,...
    'EdgeColor','none','HorizontalAlignment','left','Interpreter',...
    'none','Color',[0.1 0.1 0.1],'FontSize',8,'FitBoxToText','on');
set(h.textbox,'Position',[0.01 0.007 0.4 0.022]);
h.Ann = addExptNameToFig(h.fig,expt);
% Make figure visible
set(h.fig,'Visible','on');
if closeFig
    close(h.fig);
end
end

function setRasterPSTHpos(h)

nstim = length(h.r.ax);
ncol = ceil(nstim/2);
rrelsize = 0.65;                      % Relative size PSTH to raster
prelsize = 1-rrelsize;

% Set matrix position
margins = [0.05 0.02 0.05 0.005];
matpos = [margins(1) 1-margins(2) 0.37 1-margins(4)];  % Normalized [left right bottom top]

% Set space between plots
s1 = 0.003;
s2 = 0.035;
s3 = 0.02;

% Compute heights
rowheight = (matpos(4) - matpos(3))/2;
pheight = (rowheight-s1-s2)*prelsize;
rheight = (rowheight-s1-s2)*rrelsize;

% Compute width
width = (matpos(2)-matpos(1)-(ncol-1)*s3)/ncol;

% Row positions
p1bottom = matpos(3) + rowheight;
p2bottom = matpos(3);
r1bottom = p1bottom + pheight + s1;
r2bottom = p2bottom + pheight + s1;

% Compute complete positions
for i = 1:nstim
    if i <= ncol
        col = matpos(1)+(width+s3)*(i-1);
        p{i} = [col p1bottom width pheight];
        r{i} = [col r1bottom width rheight];
    elseif i > ncol
        col = matpos(1)+(width+s3)*(i-1-ncol);
        p{i} = [col p2bottom width pheight];
        r{i} = [col r2bottom width rheight];
    end
end

% Set positions
set([h.psth.ax; h.r.ax],'Units','normalized')
set(h.psth.ax,{'Position'},p')
set(h.r.ax,{'Position'},r')
end