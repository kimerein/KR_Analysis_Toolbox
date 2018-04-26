function [varargout]=calcMAndSForUnit_oneTrial(spikes,window,useThisTrial)

binsize=1; % in ms
% Convert binsize from ms to s
binsize = binsize/1000;
% Get count
edges=window(1):binsize:window(2);
cspikes=filtpartialspikes(spikes,0,'trials',useThisTrial);
if isempty(cspikes.spiketimes)
    n=zeros(1,length(edges)-1);
else
    n=histc(cspikes.spiketimes,edges);
end
varargout{1} = mean(n);
varargout{2}=n;
end