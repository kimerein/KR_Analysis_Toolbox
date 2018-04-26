function [ma,p1_final]=solveForP1(F_Thal,F_V1,NTC,set_ma,set_p1_final)

allF_Thal=F_Thal;
allF_V1=F_V1;

F_Thal.F1=nanmean(F_Thal.F1,1);
F_Thal.F2=nanmean(F_Thal.F2,1);
F_V1.F1=nanmean(F_V1.F1,1);
F_V1.F2=nanmean(F_V1.F2,1);

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

NTC.freqsX2=nan(size(NTC.F_NTC));
for i=1:length(NTC.freqsX2)
    findInd=find(NTC.full_freqs>=2*freqs(i),1,'first');
    NTC.freqsX2(i)=NTC.full_F_NTC(findInd);
end

% First p1 estimate
p1_1=F_V1.F1./(F_Thal.F1.*NTC.F_NTC);

% Second p1 estimate
p1_2=1-(F_V1.F2-F_Thal.F2.*NTC.freqsX2)./(F_Thal.F1.*NTC.freqsX2);

% For a estimate
a=(1./(F_Thal.F1-F_Thal.F2)).*((F_V1.F1./NTC.F_NTC)+(F_V1.F2./NTC.freqsX2));
ma=mean(a(1:11));

p1_final=(1/ma).*p1_1;

if ~isempty(set_ma)
    ma=set_ma;
end
if ~isempty(set_p1_final)
    p1_final=set_p1_final;
end
% p1_final=1-((F_V1.F2./ma)-F_Thal.F2.*NTC.freqsX2)./(F_Thal.F1.*NTC.freqsX2);

% figure();
% semilogx(freqs,p1_1,'Color','b');
% hold on; 
% semilogx(freqs,p1_2,'Color','r');
% p=mean([p1_1 p1_2]);
% 
% figure(); 
% semilogx(freqs,p*F_Thal.F1.*NTC.F_NTC,'Color',[0.5 0 0]);
% hold on; 
% semilogx(freqs,F_Thal.F1,'Color','g');
% semilogx(freqs,F_V1.F1,'Color','k');

figure(); 
semilogx(freqs,a,'Color','m');

for i=1:size(allF_Thal.F1,1)
    V1preds.F1(i,:)=ma.*p1_final.*allF_Thal.F1(i,:).*NTC.F_NTC;
    V1preds.F2(i,:)=ma.*(1-p1_final).*allF_Thal.F1(i,:).*NTC.freqsX2+ma.*allF_Thal.F2(i,:).*NTC.freqsX2;
end

plotFreqResponse(allF_Thal.F1,[],[]);
plotFreqResponse(V1preds.F1,[],[]);
plotFreqResponse(allF_V1.F1,[],[]);

plotFreqResponse(allF_Thal.F2,[],[]);
plotFreqResponse(V1preds.F2,[],[]);
plotFreqResponse(allF_V1.F2,[],[]);

% figure(); 
% semilogx(freqs,(1-p).*F_Thal.F1.*NTC.freqsX2+F_Thal.F2.*NTC.freqsX2,'Color',[0.5 0 0]);
% hold on; 
% semilogx(freqs,F_Thal.F2,'Color','g');
% semilogx(freqs,F_V1.F2,'Color','k');