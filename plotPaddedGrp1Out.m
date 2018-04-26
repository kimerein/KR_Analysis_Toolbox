function plotPaddedGrp1Out(hax,bestt,finterp,sub,freqBand,c)
 
temp=nanmean(sub(finterp>=freqBand(1) & finterp<=freqBand(2),(bestt-0.5)<1.4),1);
plot(hax,bestt(bestt-0.5<1.4)-0.5,temp,'Color',c);
hold on; 
temp2=(nanmean(sub(finterp>=freqBand(1) & finterp<=freqBand(2),(bestt-0.5)>=1.4),1)-min(nanmean(sub(finterp>=freqBand(1) & finterp<=freqBand(2),(bestt-0.5)>=1.4),1)))-(-10667.*(bestt(bestt-0.5>=1.4)-0.5-min(bestt(bestt-0.5>=1.4)-0.5)));
connect=temp(end)-temp2(1);
plot(hax,bestt(bestt-0.5>=1.4)-0.5,temp2+connect,'Color',c);
