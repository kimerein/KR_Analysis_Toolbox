function pvalue=getBinomial_pValue(nullHypothesis_pOfEvent,observedEvents)
% Given a binomial probability distribution function and a given
% probability of observing an event by chance, do the observed events lead
% us to the conclusion that we are observing a random or non-random
% process?
% INPUT PARAMETERS
% nullHypothesis_pOfEvent       a double; under the null hypothesis, probability of
%                               observing event
% observedEvents                a vector of length n where n is the number
%                               of trials; each element is 1 if the event
%                               was observed for that trial and 0 if the
%                               event was not observed for that trial
% OUTPUT
% pvalue                        a p-value for the hypothesis that
%                               observedEvents can be entirely explained by
%                               chance, given nullHypothesis_pOfEvent and a
%                               binomial distribution


pOfEvent=nullHypothesis_pOfEvent;

n=10000;
observedEvents_underNullHypothesis=zeros(n,1);
for i=1:n
    fakeObs=zeros(length(observedEvents),1);
    for j=1:length(observedEvents)
        guess=rand(1,1);
        if guess<pOfEvent
            fakeObs(j)=1;
        elseif guess>pOfEvent
            fakeObs(j)=0;
        else
            j=j-1;
        end
    end
    if sum(fakeObs)>sum(observedEvents)
        observedEvents_underNullHypothesis(i)=1;
    else
        observedEvents_underNullHypothesis(i)=0;
    end
end
pvalue=sum(observedEvents_underNullHypothesis)/n;