function fracBurst=makePreAndPostISIPlot(spikes,assignToUse,noLEDcond,showFigs)

% spikes=filtspikes(spikes,0,'assigns',assignToUse,'led',noLEDcond);
spikes=filtspikes(spikes,0,'assigns',assignToUse);

temp=[];
temp1=[];
for i=1:length(noLEDcond)
    spikes=makeTempField(spikes,'led',noLEDcond(i));
    temp(i,:)=spikes.temp;
    temp1(i,:)=spikes.sweeps.temp;
    spikes.temp=[];
    spikes.sweeps.temp=[];
end
spikes.temp=sum(temp,1)>=1;
spikes.sweeps.temp=sum(temp1,1)>=1;
spikes=filtspikes(spikes,0,'temp',1);

pres=zeros(length(spikes.spiketimes-2),1);
posts=zeros(length(spikes.spiketimes-2),1);

for i=2:length(spikes.spiketimes)-1
    pres(i-1)=spikes.spiketimes(i)-spikes.spiketimes(i-1);
    posts(i-1)=spikes.spiketimes(i+1)-spikes.spiketimes(i);
end
% for i=2:length(spikes.spiketimes)-1
%     pres(i-1)=spikes.unwrapped_times(i)-spikes.unwrapped_times(i-1);
%     posts(i-1)=spikes.unwrapped_times(i+1)-spikes.unwrapped_times(i);
% end

if showFigs==1
    if length(pres)>5000
        r=randsample(length(pres),5000);
        plotpres=pres(r);
        plotposts=posts(r);
    else
        plotpres=pres;
        plotposts=posts;
    end
    figure(); 
    maloglog(plotpres*1000,plotposts*1000);
end

total=length(pres);
bursts=sum(pres*1000<=6 & posts*1000<=6);
fracBurst=bursts/total;
