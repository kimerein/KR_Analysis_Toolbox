freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

diffs=zeros(length(freqs),1);
cons=zeros(length(freqs),1);
leds=zeros(length(freqs),1);
for i=1:length(freqs)
    [f,avSpec]=makePowerSpectrum([l1(ismember(ledConds,freqs(i)),xpoints>1 & xpoints<3); l2(ismember(ledConds2,freqs(i)),xpoints>1 & xpoints<3)],3200);
    [f2,avSpec2]=makePowerSpectrum([l1(ismember(ledConds,freqs(i)+0.05),xpoints>1 & xpoints<3); l2(ismember(ledConds2,freqs(i)+0.05),xpoints>1 & xpoints<3)],3200);
    conAmp=mean(avSpec(f>freqs(i)-0.5 & f<freqs(i)+0.5));
    ledAmp=mean(avSpec2(f2>freqs(i)-0.5 & f2<freqs(i)+0.5));
    cons(i)=conAmp;
    leds(i)=ledAmp;
    diffs(i)=conAmp-ledAmp;
end

figure(); 
plot(freqs,cons,'Color','k');
hold on; 
plot(freqs,leds,'Color','r');
plot(freqs,diffs,'Color','b');