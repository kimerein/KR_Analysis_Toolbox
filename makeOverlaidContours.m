function makeOverlaidContours(xtimes,f,allC,allCV1,n,perc)

suballC=allC(:,f>=55.5 & f<=70)'.^n;
suballCV1=allCV1(:,f>=55.5 & f<=70)'.^n;
figure(); subplot(2,1,1);
imagesc(xtimes,f(f>=55.5 & f<=70),suballC);
subplot(2,1,2);
imagesc(xtimes,f(f>=55.5 & f<=70),suballCV1);

prc1=prctile(suballC(1:end),70);
prc2=prctile(suballCV1(1:end),85);

figure(); 
subplot(3,1,1);
imagesc(xtimes,f(f>=55.5 & f<=70),suballC>=prc1);
subplot(3,1,2);
imagesc(xtimes,f(f>=55.5 & f<=70),suballCV1>=prc2);
subplot(3,1,3);
imagesc(xtimes,f(f>=55.5 & f<=70),(suballC>=prc1)+(suballCV1>=prc2));


