% script for concatenating unit data across experiments

%T1spikes=filtspikes(spikes,0,'assigns',newassignsinfo.original_assigns(newassignsinfo.trode==5));

% s1=T1234spikes;
a=unique(T1234spikes.assigns); 
[b,ind]=unique(T1234spikes.trials);
% c=0;
c=length(unitByUnitStimcond); 
for i=c+1:c+length(a)
unitByUnitTrials{i}=unique(T1234spikes.trials);
unitByUnitStimcond{i}=T1234spikes.stimcond(ind);
unitByUnitLED{i}=T1234spikes.led(ind);
end