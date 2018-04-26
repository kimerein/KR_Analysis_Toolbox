function [normCurves]=normalizeFreqResponse(responseCurves)

% Keep zero, normalize to peak
normCurves=zeros(size(responseCurves));
for i=1:size(responseCurves,1)
    currPeak=max(responseCurves(i,:));
    normCurves(i,:)=responseCurves(i,:)./currPeak;
end