function currCurve=tryLotsOfNTCs(thalcells,V1cells,bigFres)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
tryNTC=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60 70 80 90 100];

c=linspace(0,0.8,length(tryNTC));
figure();
for i=1:length(tryNTC)
    x=0:0.0001:10;
    y=exp(-x./(tryNTC(i)/1000));
    fromNTC=getTimeConstantFFT(x,y,0);
    
    [~,predOut]=doFreqAnalysis(thalcells,V1cells,fromNTC,bigFres);
    close all;
    
    [~,currCurve(i,:)]=plotFreqResponse_upperLimit(predOut,[],[],50,50,50,0,bigFres);
    
    semilogx(freqs,currCurve(i,:),'Color',[c(i) c(i) c(i)]);
end