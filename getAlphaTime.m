function [alphaTime,notAlphaTime]=getAlphaTime(alphaStates,trialDuration)

makeFigure=1;

alphaTime=zeros(1,length(alphaStates));
notAlphaTime=zeros(1,length(alphaStates));
for i=1:length(alphaStates)
    curr=alphaStates{i};
    if ~isempty(curr)
        alphaTime(i)=sum(curr(:,2)-curr(:,1));
        notAlphaTime(i)=trialDuration-alphaTime(i);
    end
end

if makeFigure==1
    figure(); 
    plot(([1:length(alphaStates)]-1)*trialDuration,alphaTime);
end