function [freqs,p]=getUnitPeakResponseAmps(allspikes,useAssign)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% freqs=[1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
% freqs=[1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05 16.05 18.05 20.05 30.05 40.05 50.05 60.05];
% useLEDcond=60;
usestimcond=1:9;
stimWindow=[1 3];
baseWindow=[0 1];
showFigs=1;    

if ~isempty(useAssign)
    allspikes=filtspikes(allspikes,0,'assigns',useAssign);
end
p=zeros(1,length(freqs));
for i=1:length(freqs)
    useLEDcond=freqs(i);
    for j=1:length(useLEDcond)
        spikes=makeTempField(allspikes,'led',useLEDcond);
        temp(j,:)=spikes.temp;
        temp1(j,:)=spikes.sweeps.temp;
        spikes.temp=[];
        spikes.sweeps.temp=[];
    end
    spikes.temp=sum(temp,1)>=1;
    spikes.sweeps.temp=sum(temp1,1)>=1;
    spikes=filtspikes(spikes,0,'temp',1,'stimcond',usestimcond);
    % spikes=filtspikes(spikes,0,'led',useLEDcond,'stimcond',usestimcond);
    s=spikes.spiketimes;
    ax=[];
    [h,ax,~,~,~,x,y]=psthForCycle_noShow(spikes,20,ax,0,4);
    if showFigs==1
        if i==1
            figure(); 
            plot(x,y);
            title('PSTH');
        end
    end
    p(i)=max(y(x>stimWindow(1) & x<=stimWindow(2)));
    pSpont(i)=max(y(x>=baseWindow(1) & x<=baseWindow(2)));
end

p=p-pSpont;
if showFigs==1
    figure(); 
    plot(freqs,p);
end
        