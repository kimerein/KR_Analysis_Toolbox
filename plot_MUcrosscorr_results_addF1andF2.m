function [alldiags,diags,diagsF2]=plot_MUcrosscorr_results_addF1andF2(MUcrosscorr)

freqs={'1';'2';'4';'6';'8';'10';'12';'14';'16';'18';'20';'30';'40';'50';'60'};

pStim=MUcrosscorr.pStim{1};

diags=zeros(1,size(pStim,1));
diagsF2=zeros(1,size(pStim,1));
for i=1:length(diags)
    diags(i)=pStim(i,i);
end
diagsF2(1)=pStim(1,2);
diagsF2(2)=pStim(2,3);
diagsF2(3)=pStim(3,5);
diagsF2(4)=pStim(4,7);
diagsF2(5)=pStim(5,9);
diagsF2(6)=pStim(6,11);
diagsF2(7)=nan;
diagsF2(8)=nan;
diagsF2(9)=nan;
diagsF2(10)=nan;
diagsF2(11)=pStim(11,13);
diagsF2(12)=pStim(12,15);
diagsF2(13)=nan;
diagsF2(14)=nan;
diagsF2(15)=nan;

diagsF2(isnan(diagsF2))=0;
alldiags=diags+diagsF2;

figure();
plot(MUcrosscorr.freqs{1},alldiags);
xlabel('Stim Freqs (Hz)');
ylabel('Response');

