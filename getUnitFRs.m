function [ledR,noLedR,sig]=getUnitFRs(spikes,useTheseAssigns)

window=[1.26 1.32];
ledOff=[0];
ledOn=[5];
stimcond=9;

ledR=zeros(length(useTheseAssigns),1);
noLedR=zeros(length(useTheseAssigns),1);
sig=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    useTheseSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',stimcond);
    [noLedR(i),temp,noledSig]=calcMeanAndStdDuringWindow(filtspikes(useTheseSpikes,0,'led',ledOff),window);
    [ledR(i),temp,ledSig]=calcMeanAndStdDuringWindow(filtspikes(useTheseSpikes,0,'led',ledOn),window);
    sig(i)=mattest(noledSig',ledSig');
end