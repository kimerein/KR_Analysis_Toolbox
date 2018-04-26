function trialsForConds=matchManyNumberOfTrials(trialsForConds)
% trialsForConds should be a cell array with n elements, where there are
% n distinct conditions
% each element of trialsForConds should be a vector specifying 
% the trial numbers for that condition

minLength=length(trialsForConds{1});
for i=1:length(trialsForConds)
    if length(trialsForConds{i})<minLength
        minLength=length(trialsForConds{i});
    end
end

for i=1:length(trialsForConds)
    theseTrials=trialsForConds{i};
    takeInds=randperm(length(theseTrials));
    trialsForConds{i}=theseTrials(sort(takeInds(1:minLength)));
end