function [data, FsDownSamp] = downSamp(y,Fs,DownSampFactor)
% Down-samples data y
%
%   Created: 2/10 - SRO
%   Modified: 4/5/10 - SRO

if size(y,2) > 1
    type = 'matrix';
else
    type = 'vector';
end

switch type
    case 'vector'
        sizeData = length(y);
        FsDownSamp = Fs/DownSampFactor;                              % Adjusted sampling frequency for filter
        data = y(1:DownSampFactor:sizeData);
    case 'matrix'
        sizeData = size(y,2);
        FsDownSamp = Fs/DownSampFactor;
        data = y(:,1:DownSampFactor:sizeData);
end



