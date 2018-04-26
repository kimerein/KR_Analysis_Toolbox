function basicFiringRate(expt,spikes,unitTag,fileInd,b,saveTag)
%
% INPUT
%   expt: Experiment struct
%   unitTag: Tag of the form 'trode_assign', e.g 'T2_15'
%   fileInd: Vector of file indices to be included in analysis.
%   b: Flag structure with field b.save, b.print, b.pause, b.close

% Created: 11/3/11 - KR
disp('Making basic firing rate figure');

figTypeName='basicFR';

% Choose binsize for PSTHs
binsize=50; % in ms

if nargin < 4
    b.pause = 0;
    b.save = 0;
    b.print = 0;
    b.close = 0;
    saveTag = '';
end

% Rig defaults
RigDef = RigDefs;

% Set cond struct
if isempty(expt.analysis.other.cond.values)
    % have not yet set expt.analysis.contrast.cond
    disp('Remember to save expt struct after entering analysis params and before analyzing!');
    return
end
cond = expt.analysis.other.cond;
if strcmp(cond.type,'no')
    spikes.all=ones(size(spikes.spiketimes));
    cond.values={1};
    cond.type='all';
else
    cond.type='led';
end

% Temporary color
if ~isfield(cond,'color')
    cond.color=[];
end
if isempty(cond.color)
    cond.color = {[0.1 0.1 0.1],[1 0.25 0.25],[0 0 1],[1 0.5 0],[1 0 1],[0.3 0.3 0.3],[0.7 0.7 0.7]};
end
gray = [0.6 0.6 0.6];
green = [0 1 0];
blue = [0 120/255 200/255];
red = [1 0.25 0.25];

% Set time window struct
if ~isfield(expt.analysis.other,'w')
    expt.analysis.other.w=[];
end
if ~isfield(expt.analysis.other.w,'spont')
    % Prompt user for window information
    prompt={'Baseline window:','Stimulus window:','Onset response window:','Offset response window:','LED On window:'};
    dlg_title='Get special time windows';
    num_lines = 1;
    %def = {'0 1','1 4','1 1.5','4 4.5','1 2'};
    def = {'0 5','5 5','0 0.5','4.5 5','1 1.25'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    expt.analysis.other.w.spont=str2num(answer{1});
    expt.analysis.other.w.stim=str2num(answer{2});
    expt.analysis.other.w.on=str2num(answer{3});
    expt.analysis.other.w.off=str2num(answer{4});
    expt.analysis.other.w.ledOn=str2num(answer{5});
    w=expt.analysis.other.w;
else
    w=expt.analysis.other.w;
end

% Figure layout
h.fig = linearArrayFigSetup;
set(h.fig,'Visible','off');

% Set save name suffix (saveTag)
if isempty(saveTag)
    saveTag = [unitTag '_basicFR'];
else
    saveTag = [unitTag '_' saveTag];
end

% Set expt struct as appdata
setappdata(h.fig,'expt',expt);
setappdata(h.fig,'figText',saveTag);

% Get trial duration
getTrialDurationFromUser=0;
duration=0;
if isstruct(spikes.info.detect)
    duration=spikes.info.detect.dur(1);
else
    duration=expt.files.duration(fileInd(1));
end

% Add save figure button
addSaveFigTool(h.fig);

% Extract spikes for files
spikes = filtspikes(spikes,0,'fileInd',fileInd);

if isempty(spikes.spiketimes) % Check that there are spikes
    return
end

% Needs to be earlier, so event_channel is filtered, because
% spikes.info.detect.event_channel is not filtered by filterspikes
% if length(spikes.info.detect)==1 % Not concatenated spikes
%     spikes.event_channel=spikes.info.detect.event_channel;
%     spikes.event_channel=spikes.event_channel';
% end

% Set NaNs = 0
spikes.led(isnan(spikes.led)) = 0;
spikes.sweeps.led(isnan(spikes.sweeps.led)) = 0;
if strcmp(cond.type,'led')
    for j=1:length(cond.values)
        tempSpikes=makeTempField(spikes,'led',cond.values{j});
        cspikes(1,j)=filtspikes(tempSpikes,0,'temp',1);
    end
else
    cspikes(1,1)=spikes;
end

% Make PSTH for each cspikes substruct
for m=1:size(cspikes,1)
    h.psth.ax(m)=axes;
    for n=1:size(cspikes,2)
        [h.psth.lines(m,n) temp h.psth.n(m,:,n) centers] = psth(cspikes(m,n),binsize,h.psth.ax(m),1,duration);
    end
    %setTitle(h.psth.ax(m),strcat('Ch',num2str(eventChs(m))),7);
    %set(h.psth.ax,'Title',strcat('Ch',num2str(eventChs(m))),'Position',[1 0 1]);
    %set(get(h.psth.ax(m),'Title'),'String',strcat('Ch',num2str(eventChs(m))),'Position',[0.4 7 1]);
end
h.psth.ax=h.psth.ax';

% Add handles to appropriate condition field
for n = 1:size(cspikes,2)
    h.(cond.tags{n}) = [];
end
for n = 1:size(cspikes,2)
    h.(cond.tags{n}) = [h.(cond.tags{n}); h.psth.lines(:,n)];
end


% Set axes properties
ymax=setSameYmax(h.psth.ax,15);
for i=1:length(h.psth.ax)
    addStimulusBar(h.psth.ax(i),[w.stim ymax]);
    addStimulusBar(h.psth.ax(i),[w.ledOn ymax*0.97],'','red');
end
removeAxesLabels(h.psth.ax(1:length(h.psth.ax)-1));
defaultAxes(h.psth.ax,0.05,0.05);
ylabel('spikes/s');

% Compute average firing rate for different windows
% for i=1:size(cspikes,1)
%     for j=1:size(cspikes,2)
%         wnames=fieldnames(w);
%         for n=1:length(wnames)
%             temp=wnames{n};
%             [fr.(temp)(i,j,1) fr.(temp)(i,j,2)]=computeFR(cspikes(i,j),w.(temp)); % average and SEM
%         end
%     end
% end
% 
% % Make category plot for each time window
% for j=1:size(cspikes,1)
%     for i = 1:length(wnames)
%         temp = wnames{i};
%         [h.avgfr.(temp).l(j) h.avgfr.(temp).ax(j)] = plotCategories(fr.(temp)(j,:,1),cond.tags,fr.(temp)(j,:,2),'');
%         setTitle(gca,temp,7);
%     end
% end
% defaultAxes(h.avgfr.spont.ax(1),0.1,0.48);

% Set axes locations
setLinearArrayPositions(h);

% Set colors
for i = 1:length(cond.tags)
    set(h.(cond.tags{i}),'Color',cond.color{i})
end

exptInfo = unitTag;
h.textbox = annotation('textbox',[0 0 0.3 0.022],'String',exptInfo,...
    'EdgeColor','none','HorizontalAlignment','left','Interpreter',...
    'none','Color',[0.1 0.1 0.1],'FontSize',8,'FitBoxToText','on');
set(h.textbox,'Position',[0.01 0.013 0.4 0.022]);

set(h.fig,'Visible','on')

if b.pause
    reply = input('Do you want to print? y/n [n]: ', 's');
    if isempty(reply)
        reply = 'n';
    end
    if strcmp(reply,'y')
        b.print = 1;
    end
end
sname = [RigDef.Dir.Fig expt.name '_' unitTag '_' figTypeName];
if b.save
    disp(['Saving' ' ' sname])
    saveas(h.fig,sname,'fig')
    saveas(h.fig,sname,'epsc')
    
%     format long
%     format compact
    set(0, 'defaultFigurePaperType', 'A4')
    set(0, 'defaultFigurePaperUnits', 'centimeters')
    set(0, 'defaultFigurePaperPositionMode', 'auto')
    saveas(h.fig,sname,'pdf')
    
    sname = [sname '.epsc'];
    export_fig sname
end
if b.print
    print('-dwinc',h.fig)
    disp(['Printing' ' ' sname])
end
if b.close
    close(h.fig)
end
end

function setLinearArrayPositions(h)

h.mat(1).params.matpos=[0.02 0.02 0.94 0.93]; % [left top width height]
h.mat(1).params.figmargin=[0 0 0 0.05];
h.mat(1).params.matmargin=[0 0 0 0];
h.mat(1).params.cellmargin=[0.05 0.05 0.005 0.005];
h.mat(1).ncol=1;
h.mat(1).nrow=length(h.psth.ax);
for i=1:length(h.psth.ax)
    h.mat(1).h(i)=h.psth.ax(i);
end

% h.mat(2).params.matpos=[0.61 0 0.385 1];
% h.mat(2).params.figmargin=[0 0 0 0];
% h.mat(2).params.matmargin=[0 0 0 0];
% h.mat(2).params.cellmargin=[0.03 0.03 0.03 0.03];
% h.mat(2).ncol=5;
% h.mat(2).nrow=length(h.psth.ax);
% for i=1:length(h.psth.ax)
%     h.mat(2).h(((i-1)*5)+1)=h.avgfr.spont.ax(i);
%     h.mat(2).h(((i-1)*5)+2)=h.avgfr.stim.ax(i);
%     h.mat(2).h(((i-1)*5)+3)=h.avgfr.on.ax(i);
%     h.mat(2).h(((i-1)*5)+4)=h.avgfr.off.ax(i);
%     h.mat(2).h(((i-1)*5)+5)=h.avgfr.ledOn.ax(i);
% end

% Place axes on axesmatrix
for i = 1:length(h.mat)
    ind = 1:length(h.mat(i).h);
    setaxesOnaxesmatrix(h.mat(i).h,h.mat(i).nrow,h.mat(i).ncol,ind, ...
        h.mat(i).params,h.fig);
end

end

function hFig=linearArrayFigSetup(varargin)

if nargin>0
    hFig=varargin{1};
else
    hFig=[];
end

if ~isempty(hFig)
    hFig=figure(hFig,'Visible','off');
else
    hFig=figure('Visible','off');
end

% Set monitor position of figure, if you want
% Set paper position in inches (where figures will be printed on paper)
leftMargin=0.05;
rightMargin=0.05;
bottomMargin=0.1;
topMargin=0.05;

% Monitor position
% monitorPosition=[20 0 500 880];
%monitorPosition=[20 0 500 500];
monitorPosition=[20 0 500 500];

PaperPosition=[leftMargin bottomMargin (11-leftMargin-rightMargin) (8.5-topMargin-bottomMargin)];
set(hFig,'PaperPositionMode','auto','PaperOrientation','portrait','PaperPosition',PaperPosition,'Color',[0.94 0.94 0.94],'Position',monitorPosition);

% Make figure visible
set(hFig,'Visible','on');

end

        