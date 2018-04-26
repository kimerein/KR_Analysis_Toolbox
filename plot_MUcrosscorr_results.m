function [diags,pStimMinusOffDiags]=plot_MUcrosscorr_results(MUcrosscorr)

subtractOffDiags=0;
transformToAmps=1;

freqs={'1';'2';'4';'6';'8';'10';'12';'14';'16';'18';'20';'30';'40';'50';'60'};

% figure();
% imagesc(MUcrosscorr.pStim{1});
% set(gca,'XTick',[1:size(MUcrosscorr.pStim{1},1)]);
% set(gca,'XTickLabel',freqs);
% set(gca,'YTick',[1:size(MUcrosscorr.pStim{1},1)]);
% set(gca,'YTickLabel',freqs);
% xlabel('Response Freqs (Hz)');
% ylabel('Stim Freqs (Hz)');
% title('MU Cross-Corr During Stim');

% % figure();
% % imagesc(MUcrosscorr.pSpont{1});
% % set(gca,'XTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'XTickLabel',freqs);
% % set(gca,'YTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'YTickLabel',freqs);
% % xlabel('Response Freqs (Hz)');
% % ylabel('Stim Freqs (Hz)');
% % title('MU Cross-Corr During Spont');
% 
% % figure();
% % imagesc(MUcrosscorr.pStim{1}-MUcrosscorr.pSpont{1});
% % set(gca,'XTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'XTickLabel',freqs);
% % set(gca,'YTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'YTickLabel',freqs);
% % xlabel('Response Freqs (Hz)');
% % ylabel('Stim Freqs (Hz)');
% % title('MU Cross-Corr Stim Minus Spont');

if transformToAmps==1
    pStim=sqrt(MUcrosscorr.pStim{1});
else
    pStim=MUcrosscorr.pStim{1};
end
offDiags=zeros(1,size(pStim,2));
for i=1:size(pStim,2)
    currRow=i;
%     offDiags(i)=mean(pStim([1:currRow-1 currRow+1:end],i));
    offDiags(i)=mean(pStim(i,[1:currRow-1 currRow+1:end]));
end
offDiagMatrix=repmat(offDiags,size(pStim,1),1);
pStimMinusOffDiags=pStim-offDiagMatrix;
% % figure();
% % imagesc(pStimMinusOffDiags);
% % set(gca,'XTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'XTickLabel',freqs);
% % set(gca,'YTick',[1:size(MUcrosscorr.pStim{1},1)]);
% % set(gca,'YTickLabel',freqs);
% % xlabel('Response Freqs (Hz)');
% % ylabel('Stim Freqs (Hz)');
% % title('MU Cross-Corr Stim Minus Off-Diagonals');
if subtractOffDiags==0
    pStimMinusOffDiags=pStim;
end

diags=zeros(1,size(pStimMinusOffDiags,1));
for i=1:length(diags)
    diags(i)=pStimMinusOffDiags(i,i);
%     diags(i)=pStim(i,i);
end
% figure();
% semilogx(MUcrosscorr.freqs{1},diags);
% xlabel('Stim Freqs (Hz)');
% ylabel('Response');

