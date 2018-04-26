function [hPsth,maxVal,binWidth,ONresponse,OFFresponse]=K_superPSTH(spikes,idealOrSet,binWidth,colorMap,colorOffset,params,tit,ax,quantResponse,ONresponseWindow,OFFresponseWindow) 
% Plots PSTH for different stimulus conditions               
% EACH STIMULUS CONDITION IS SUPERIMPOSED ON SAME PLOT
% 
% PARMATERS:
% spikes is a vector of spikes structures for different stimulus conditions
% each element of spikes contains a structure for a different stimulus
% condition
%
% idealOrSet is method for choosing PSTH bin width ...
% idealOrSet is 'ideal' -> bin width is that which gives the maximal
%   range of firing rates in the PSTH
% idealOrSet is 'set' -> bin width is set by the user in variable binWidth
% 
% colorMap is the color code for different stimulus conditions
% each three-element row of colorMap is the color for a different stimulus
% condition
% PSTH lines will be colored according to stimulus condition
% Alternatively, pass in a three-element vector specifying color for all
% lines, e.g., colorMap=[0 0 0]
%
% params are the parameters for this stimulus set (i.e., stim. ON period 
% duration, baseline measure duration, total trial length, etc.)
% the params structure must AT LEAST CONTAIN:
% params.ONstart = the time the stimulus appeared on the screen
% params.totalTrialLength = the total acquired time for each sweep
% params.ONlength = the duration of the stimulus
%
% tit is the string to display at the top of the PSTH
% 
% ax is a four-element specifying [x_axis_lower_limit x_axis_upper_limit
%                                  y_axis_lower_limit y_axis_upper_limit]
% if ax is the empty array, the PSTH axes are automatically chosen
% 
% Returns the ON and OFF responses quantified by quantMethod
%
% Time windows to use for the calculation of ON and OFF responses are given
% by ONresponseWindow and OFFresponseWindow, respectively
%
% ONresponseWindow is a two-element vector specifying the beginning and end
% of the time period to consider for ON response; same for
% OFFresponseWindow
% 
% quantResponse is the method for quantifying the ON or OFF response
%       can be 'peak','integral' or 'average'
%       Method - Response is the ... 
%       'peak' - peak firing rate during the specified period
%       'integral' - area under the firing rate curve of the specified
%       period
%       'average' - average of the firing rate across the entire specified
%       period
%
% RETURNS:
% returns handle to figure and maximum firing rate across all stimulus
% conditions; also, binWidth and calculated ON and OFF responses

% Get current figure and current axes
hax = gca;

% Default hPsth
hPsth=hax;

% Default response values
ONresponse=[];
OFFresponse=[];

if strcmp(idealOrSet,'set')
    % Use user-provided bin width
    % binWidth is in s
    b=binWidth;
    edges=0:b:params.totalTrialLength;
    x=zeros(1,length(edges));
    maxVal=0.1;
    for i=1:length(spikes)   
        if isempty(spikes(i).spiketimes)
            ps=zeros(1,length(x));
            continue
        end
        trials=spikes(i).sweeps.trials;
        ps=histc(spikes(i).spiketimes,edges);
        ps=ps/(b*length(trials));
        for j=1:length(edges)-1
            x(j)=mean([edges(j) edges(j+1)]);
        end
        x(end)=mean([edges(end) params.totalTrialLength]);
        if size(colorMap,1)==1
            %hPsth=line(edges-params.ONstart,ps,'Parent',hax,'LineWidth',1,'Color',colorMap);
            hPsth=line(edges,ps,'Parent',hax,'LineWidth',1.5,'Color',colorMap);
        else
            if colorOffset==0
                c=colorMap(spikes(i).stimcond(1),:);
            else
                c=colorMap(colorOffset+i,:);
            end
            %hPsth=line(edges-params.ONstart,ps,'Parent',hax,'LineWidth',1.5,'Color',c);  
            hPsth=line(edges,ps,'Parent',hax,'LineWidth',1.5,'Color',c);  
        end
        hold on;
        m=max(ps);
        if m>maxVal
            maxVal=m;
        end
        ONresponse=[ONresponse; max(ps(intersect(find(edges>ONresponseWindow(1)),find(edges<ONresponseWindow(2)))))];
        OFFresponse=[OFFresponse; max(ps(intersect(find(edges>OFFresponseWindow(1)),find(edges<OFFresponseWindow(2)))))];
    end
    if ~isempty(ax)
        axis(hax,ax);
    else
        axis(hax,[0 params.totalTrialLength 0 maxVal]);
    end
    % Quantify ON and OFF responses
    switch quantResponse
        case 'peak'
%             psth_inds=find(x>ONresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
%             ONresponse=max(ps(psth_inds(i1)));
%             
%             psth_inds=find(x>OFFresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
%             OFFresponse=max(ps(psth_inds(i1)));
%             if OFFresponse>0
%                 'hi'
%             end
        case 'integral'
              for i=1:length(spikes)
                  %ONresponse=[ONresponse; getResponse(spikes(i),ONresponseWindow)];
                  %ONresponse=[ONresponse; getResponse(spikes(i),ONresponseWindow)];
              end
%             psth_inds=find(x>ONresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
%             ONresponse=sum(ps(psth_inds(i1)));
              for i=1:length(spikes)
                  %OFFresponse=[OFFresponse; getResponse(spikes(i),OFFresponseWindow)];
                  %OFFresponse=[OFFresponse; getResponse(spikes(i),OFFresponseWindow)];
              end
%             psth_inds=find(x>OFFresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
%             OFFresponse=sum(ps(psth_inds(i1)));
        case 'average'
%             psth_inds=find(x>ONresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
%             ONresponse=mean(ps(psth_inds(i1)));
%             
%             psth_inds=find(x>OFFresponseWindow(1));
%             [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
%             OFFresponse=mean(ps(psth_inds(i1)));
    end      
    %title(hax,{tit; ['Stimulus is ' num2str(params.ONlength) ' seconds in duration']});
    title(hax,tit);
%     xlabel(hax,'Time from Stimulus Onset (s)');
%     ylabel(hax,'F.R.(Hz)');
else
    % Find ideal bin width for plotting PSTH
    % "Ideal" = binning that maximizes the difference between the PSTH
    % minimum and the PSTH maximum
    % Try several bin sizes
    %binSizes=[0.05 0.1 0.2 0.3 0.5];
    binSizes=[0.1 0.2 0.3 0.5];
    maxDiff=zeros(1,5);
    x=0; % default value for x
    for i=1:length(binSizes)
        trials=spikes(1).sweeps.trials;
        edges=0:i:params.totalTrialLength;
        ps=histc(spikes(1).spiketimes,edges);
        ps=ps/(i*length(trials));
        maxDiff(i)=max(ps)-min(ps);
    end
    [m,mInd]=max(maxDiff);
    b=binSizes(mInd);
    maxVal=0.1;
    for i=1:length(spikes)
        if isempty(spikes(i).spiketimes)
            ps=0;
            continue
        end
        trials=spikes(i).sweeps.trials;
        edges=0:b:params.totalTrialLength;
        ps=histc(spikes(i).spiketimes,edges);
        ps=ps/(b*length(trials));
        x=zeros(1,length(ps));
        for j=1:length(edges)-1
            x(j)=mean([edges(j) edges(j+1)]);
        end
        x(end)=mean([edges(end) params.totalTrialLength]);
        if size(colorMap,1)==1
            %hPsth=line(edges-params.ONstart,ps,'Parent',hax,'LineWidth',1,'Color',colorMap);
            hPsth=line(edges,ps,'Parent',hax,'LineWidth',1.5,'Color',colorMap);
        else
            if colorOffset==0
                c=colorMap(spikes(i).stimcond(1),:);
            else
                c=colorMap(colorOffset+i,:);
            end
            hPsth=line(edges,ps,'Parent',hax,'LineWidth',1.5,'Color',c);
        end
        hold on;
        m=max(ps);
        if m>maxVal
            maxVal=m;
        end
        ONresponse=[ONresponse; max(ps(intersect(find(edges>ONresponseWindow(1)),find(edges<ONresponseWindow(2)))))];
        OFFresponse=[OFFresponse; max(ps(intersect(find(edges>OFFresponseWindow(1)),find(edges<OFFresponseWindow(2)))))];
    end
    if ~isempty(ax)
        axis(hax,ax);
    else
        axis(hax,[0 params.totalTrialLength 0 maxVal]);
    end
    % Quantify ON and OFF responses
    switch quantResponse
        case 'peak'
            % Need to fill in
        case 'integral'
            for i=1:length(spikes)
                %ONresponse=[ONresponse; getResponse(spikes(i),ONresponseWindow)];
            end
            for i=1:length(spikes)
                %OFFresponse=[OFFresponse; getResponse(spikes(i),OFFresponseWindow)];
            end
        case 'average'
            % Need to fill in
    end  
    %title(hax,{tit; ['Stimulus is ' num2str(params.ONlength) ' seconds in duration']});
    title(hax,tit);
%     xlabel(hax,'Time from Stimulus Onset (s)');
%     ylabel(hax,'F.R.(Hz)');
end
if isempty(ONresponse)
    ONresponse=0;
end
if isempty(OFFresponse)
    OFFresponse=0;
end