function paired=matchUnits(units1,units2,base1,base2,spikes1,spikes2,origUnitInfo1,origUnitInfo2)

maxUnitAssigns=200;

if ~isempty(base1)
    units1.led1Spikes=units1.led1Spikes-base1.led1Spikes;
    units1.led2Spikes=units1.led2Spikes-base1.led2Spikes;
end
if ~isempty(base2)
    units2.led1Spikes=units2.led1Spikes-base2.led1Spikes;
    units2.led2Spikes=units2.led2Spikes-base2.led2Spikes;
end
if ~isfield(units1,'assigns')
    if ~isempty(origUnitInfo1)
    elseif isempty(spikes1)
        disp('Error: units1 has no assigns');
        return
    else
        units1.assigns=unique(spikes1.assigns);
    end
end
if ~isfield(units2,'assigns')
    if ~isempty(origUnitInfo2)
    elseif isempty(spikes2)
        disp('Error: units2 has no assigns');
        return
    else
        units2.assigns=unique(spikes2.assigns);
    end
end
if ~isempty(origUnitInfo1)
    units1.assigns=[];
    ts=unique(origUnitInfo1.trode);
    for i=1:length(ts)
        currTassigns=origUnitInfo1.original_assigns(origUnitInfo1.trode==ts(i));
        sortedTassigns=sort(currTassigns);
        units1.assigns=[units1.assigns sortedTassigns+((i-1)*maxUnitAssigns)];
    end
end
if ~isempty(origUnitInfo2)
    units2.assigns=[];
    ts=unique(origUnitInfo2.trode);
    for i=1:length(ts)
        currTassigns=origUnitInfo2.original_assigns(origUnitInfo2.trode==ts(i));
        sortedTassigns=sort(currTassigns);
        units2.assigns=[units2.assigns sortedTassigns+((i-1)*maxUnitAssigns)];
    end
end

units1_assigns=units1.assigns;
units2_assigns=units2.assigns;
[maLength,ind]=max([length(units1_assigns) length(units2_assigns)]);

if ind==2
    temp=units1;
    units1=units2;
    units2=temp;
    switched=1;
else
    switched=0;
end
j1=1;
j2=1;
paired.led1Spikes=zeros(maLength,2);
paired.led2Spikes=zeros(maLength,2);
paired.assigns=zeros(maLength,1);
for i=1:maLength
  if j1>length(units1.assigns)
      currAss1=-1;
  else
      currAss1=units1.assigns(j1);
  end
  if j2>length(units2.assigns)
      currAss2=-10;
  else
      currAss2=units2.assigns(j2);
  end
  if currAss1==currAss2
      paired.led1Spikes(i,1)=units1.led1Spikes(j1);
      paired.led1Spikes(i,2)=units2.led1Spikes(j2);
      paired.led2Spikes(i,1)=units1.led2Spikes(j1);
      paired.led2Spikes(i,2)=units2.led2Spikes(j2);
      paired.assigns(i)=currAss1;
      j1=j1+1;
      j2=j2+1;
  else
      if ~ismember(currAss1,units2.assigns) && currAss1>0
          paired.led1Spikes(i,1)=units1.led1Spikes(j1);
          paired.led1Spikes(i,2)=0;
          paired.led2Spikes(i,1)=units1.led2Spikes(j1);
          paired.led2Spikes(i,2)=0;
          paired.assigns(i)=currAss1;
          j1=j1+1;
      end
      if ~ismember(currAss2,units1.assigns) && currAss2>0
          paired.led1Spikes(i,1)=0;
          paired.led1Spikes(i,2)=units2.led1Spikes(j2);
          paired.led2Spikes(i,1)=0;
          paired.led2Spikes(i,2)=units2.led2Spikes(j2);
          paired.assigns(i)=currAss2;
          j2=j2+1;
      end 
  end
end

if switched==1
    figure();
    scatter(paired.led1Spikes(:,2),paired.led1Spikes(:,1));
    title('Led Cond 1');
    figure(); 
    scatter(paired.led2Spikes(:,2),paired.led2Spikes(:,1));
    title('Led Cond 2');
else
    figure();
    scatter(paired.led1Spikes(:,1),paired.led1Spikes(:,2));
    title('Led Cond 1');
    figure(); 
    scatter(paired.led2Spikes(:,1),paired.led2Spikes(:,2));
    title('Led Cond 2');
end
      