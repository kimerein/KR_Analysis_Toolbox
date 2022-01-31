function normedOut_all=combineTheseNormedToQuiescentControl(filenames,whichInFilename,whichForNormVals)

normedOut_all=[];
for i=1:length(filenames)
    a=load(whichForNormVals{i});
    condenom=a.denom;
    a=load(filenames{i});
    numer=a.numer(:,whichInFilename(i));
    normedOut_part=numer./nanmean(condenom,2);
    normedOut_all=[normedOut_all; normedOut_part];
end

% for sf 0.03
% normedOut_all=[];
% for i=1:length(filenames)
%     a=load(whichForNormVals{i});
%     condenom=a.denom;
%     a=load(filenames{i});
%     numer=a.numer;
%     normedOut_part=nanmean(numer./condenom,2);
%     normedOut_all=[normedOut_all; normedOut_part];
% end