function hErr = addErrBar_asymmetric(x,y,errLow,errHigh,errType,hax,hl,col)
% function hErr = addErrBar(x,y,err,errType,hax,hl)
%
% INPUT
%   x: XData for data line
%   y: YData for data line
%   err: Error bar values
%   errType: Add error bars to either 'x' or 'y' values
%   hax: Handle to parent axes for data line
%   hl: Handle to data line
%
% OUTPUT
%   hErr: Handle to error bar line

% Created: SRO - 5/25/11
% Modified: KR - 10/16/2014


if nargin < 5 || isempty(errType)
    errType = 'y';
end

if nargin < 6 || isempty(hax)
    hax = gca;
end

if nargin < 7 || isempty(hl)
    hl = [];
end

if nargin < 8 || isempty(col)
    col = 'k';
end

if strcmp(errType,'x')
    xtmp = x;
    x = y;
    y = xtmp;
end

% Compute +/- error value
yp = y + errHigh;
ym = y - errLow;

% Make y values vector
ny = NaN(3*length(yp),1);
ny(1:3:end) = ym;
ny(2:3:end) = yp;

% Make x values vector
nx = NaN(size(ny));
nx(1:3:end) = x;
nx(2:3:end) = x;

% Plot error bars
switch  errType
    case 'y'
        hErr = line('XData',nx,'YData',ny,'Parent',hax,'Color',col);
    case 'x'
        hErr = line('XData',ny,'YData',nx,'Parent',hax,'Color',col);
end

% Format error bars
if ~isempty(hl)
    set(hErr,'Color',get(hl,'Color'));
end
