function [h,maxVal,binWidth,ONresponse,OFFresponse]=superPSTH(hax, spikes, idealOrSet, binWidth, colorMap, params, tit, ax, quantResponse, ONresponseWindow, OFFresponseWindow) 
% Plots PSTH for different stimulus conditions               
% EACH STIMULUS CONDITION IS SUPERIMPOSED ON SAME PLOT
% 
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
% returns handle to figure and maximum firing rate across all stimulus
% conditions


h=gcf;
if strcmp(idealOrSet,'set')
    % Use user-provided bin width
    b=binWidth;
    maxVal=0.1;
    for i=1:length(spikes)        
        [hi,x]=hist(spikes(i).spiketimes,b/2:b:params.totalTrialLength);
        psth=hi/b;
        plot(hax,x-params.ONstart,psth,'Color',colorMap(i,:));
        hold on;
        m=max(psth);
        if m>maxVal
            maxVal=m;
        end
    end
    if ~isempty(ax)
        axis(hax,ax);
    else
        axis(hax,[0 params.totalTrialLength 0 maxVal]);
    end
    % Quantify ON and OFF responses
    switch quantResponse
        case 'peak'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=max(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=max(psth(psth_inds(i1)));
        case 'integral'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=sum(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=sum(psth(psth_inds(i1)));
        case 'average'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=mean(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=mean(psth(psth_inds(i1)));
    end      
    title(hax,{tit; ['Stimulus is ' num2str(params.ONlength) ' seconds in duration']});
    xlabel(hax,'Time from Stimulus Onset (s)');
    ylabel(hax,'F.R.(Hz)');
else
    % Find ideal bin width for plotting PSTH
    % "Ideal" = binning that maximizes the difference between the PSTH
    % minimum and the PSTH maximum
    % Try several bin sizes
    binSizes=[0.05 0.1 0.2 0.3 0.5];
    maxDiff=zeros(1,5);
    for i=1:length(binSizes)
        b=binSizes(i);
        length(spikes)
        [psth,x]=hist(spikes(1).spiketimes,b/2:b:params.totalTrialLength);
        maxDiff(i)=max(psth)-min(psth);
    end
    [m,mInd]=max(maxDiff);
    b=binSizes(mInd);
    maxVal=0.1;
    for i=1:length(spikes)
        [hi,x]=hist(spikes(i).spiketimes,b/2:b:params.totalTrialLength);
        psth=hi/b;
        plot(hax,x-params.ONstart,psth,'Color',colorMap(spikes(i).stimcond(1),:));
        hold on;
        m=max(psth);
        if m>maxVal
            maxVal=m;
        end
    end
    if ~isempty(ax)
        axis(hax,ax);
    else
        axis(hax,[0 params.totalTrialLength 0 maxVal]);
    end
    % Quantify ON and OFF responses
    switch quantResponse
        case 'peak'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=max(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=max(psth(psth_inds(i1)));
        case 'integral'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=sum(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=sum(psth(psth_inds(i1)));
        case 'average'
            psth_inds=find(x>ONresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>ONresponseWindow(1)),x(x<ONresponseWindow(2)));
            ONresponse=mean(psth(psth_inds(i1)));
            
            psth_inds=find(x>OFFresponseWindow(1));
            [x_intersect,i1,i2]=intersect(x(x>OFFresponseWindow(1)),x(x<OFFresponseWindow(2)));
            OFFresponse=mean(psth(psth_inds(i1)));
    end  
    title(hax,{tit; ['Stimulus is ' num2str(params.ONlength) ' seconds in duration']});
    xlabel(hax,'Time from Stimulus Onset (s)');
    ylabel(hax,'F.R.(Hz)');
end