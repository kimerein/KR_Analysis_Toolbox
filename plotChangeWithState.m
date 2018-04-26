function [spontSteady,evSteady,evsubSteady]=plotChangeWithState(spont,ev,evsub,offAt,offAtEv)

% Spont
for i=1:length(spont)
    spontSizes(i)=length(spont{i});
end
arraySecondHalf=max(spontSizes-offAt);
arraySize=max(offAt)+arraySecondHalf;
spontShifted=zeros(length(spont),arraySize);
for i=1:length(spont)
    spontShifted(i,:)=[nan(1,max(offAt)-offAt(i)) spont{i} nan(1,arraySize-(length(spont{i})+max(offAt)-offAt(i)))];
end
figure(); 
plot(spontShifted');
hold all; 
plot(mean(spontShifted,1)','k');
for i=1:length(spont)
    curr=spontShifted(i,1:max(offAt));
    spontSteady(i,1)=mean(curr(~isnan(curr)));  
    curr=spontShifted(i,max(offAt)+1:end);
    spontSteady(i,2)=mean(curr(~isnan(curr)));
end

% Ev
for i=1:length(ev)
    evSizes(i)=length(ev{i});
end
arraySecondHalf=max(evSizes-offAtEv);
arraySize=max(offAtEv)+arraySecondHalf;
evShifted=zeros(length(ev),arraySize);
for i=1:length(ev)
    evShifted(i,:)=[nan(1,max(offAtEv)-offAtEv(i)) ev{i} nan(1,arraySize-(length(ev{i})+max(offAtEv)-offAtEv(i)))];
end
figure(); 
plot(evShifted');
hold all; 
plot(mean(evShifted,1)','k');
for i=1:length(ev)
    curr=evShifted(i,1:max(offAtEv));
    evSteady(i,1)=mean(curr(~isnan(curr)));
    curr=evShifted(i,max(offAtEv)+1:end);
    evSteady(i,2)=mean(curr(~isnan(curr)));
end

% EvSub
for i=1:length(evsub)
    evsubSizes(i)=length(evsub{i});
end
arraySecondHalf=max(evsubSizes-offAtEv);
arraySize=max(offAtEv)+arraySecondHalf;
evsubShifted=zeros(length(evsub),arraySize);
for i=1:length(evsub)
    evsubShifted(i,:)=[nan(1,max(offAtEv)-offAtEv(i)) evsub{i} nan(1,arraySize-(length(evsub{i})+max(offAtEv)-offAtEv(i)))];
end
figure(); 
plot(evsubShifted');
hold all; 
plot(mean(evsubShifted,1)','k');
for i=1:length(evsub)
    curr=evsubShifted(i,1:max(offAtEv));
    evsubSteady(i,1)=mean(curr(~isnan(curr)));
    curr=evsubShifted(i,max(offAtEv)+1:end);
    evsubSteady(i,2)=mean(curr(~isnan(curr)));
end
