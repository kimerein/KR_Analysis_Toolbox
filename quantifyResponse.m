function [ONresponse, OFFresponse]=quantifyResponse(spikes, )
% Returns the ON and OFF responses quantified by quantMethod
% Time windows to use for the calculation of ON and OFF responses are given
% by ONresponseWindow and OFFresponseWindow, respectively
% ONresponseWindow is a two-element vector specifying the beginning and end
% of the time period to consider for ON response; same for
% OFFresponseWindow
% quantResponse is the method for quantifying the ON or OFF response
%       can be 'peak','integral' or 'average'
%       Method - Response is the ... 
%       'peak' - peak firing rate during the specified period
%       'integral' - area under the firing rate curve of the specified
%       period
%       'average' - average of the firing rate across the entire specified
%       period

% Quant ON response
binWidth=ONresponseWindow(2)-ONresponseWindow(1);
centerBin=ONresponseWindow(1)+binWidth/2;
[hi,x]=hist(spikes.spiketimes,[centerBin-binWidth centerBin centerBin+binWidth]);
ONresponse=hi(2)/binWidth;



[hi,x]=hist(spikes(i).spiketimes,b/2:b:params.totalTrialLength);
        psth=hi/b;
        plot(x-params.ONstart,psth,'Color',colorMap(i,:));