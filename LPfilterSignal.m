function s=LPfilterSignal(signal,butterCutoff)

if isempty(butterCutoff)
%     butterCutoff=0.025;
    butterCutoff=0.05;
end

[B,A]=butter(5,butterCutoff,'low');
xFiltered=filtfilt(B,A,signal);
s=xFiltered;