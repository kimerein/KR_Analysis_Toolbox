function avFracBursts=getAverageFractionBursts(allspikespsth,burstpsth,outdir)

avFracBursts=nan(1,length(allspikespsth.psths));
for i=1:length(allspikespsth.psths)
    p=allspikespsth.psths{i};
    burstp=burstpsth.psths{i};
    allspikessum=sum(sum(p,1),2);
    burstsum=sum(sum(burstp,1),2);
    avFracBursts(i)=burstsum/allspikessum;
end

save([outdir '\avFracBursts.mat'],'avFracBursts');