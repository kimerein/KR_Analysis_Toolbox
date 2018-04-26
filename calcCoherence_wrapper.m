function [freqs,coh,coh_matrix,coh_matrix_minusOffDiags,diags]=calcCoherence_wrapper(spikes)

onlyDoDiags=1;

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% leds=[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050];
leds=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

coh_matrix_minusOffDiags=[];
diags=[];
coh=zeros(1,length(freqs));
cerr1=zeros(1,length(freqs));
cerr2=zeros(1,length(freqs));
coh_matrix=zeros(length(leds),length(freqs));
for i=1:length(freqs)
    for j=1:length(leds)
        if onlyDoDiags==1 && i~=j
            continue
        end
        [c,c1,c2]=calcCoherence(freqs(i),spikes,leds(j));
        coh_matrix(j,i)=c;
        if i==j
            coh(i)=c;
            cerr1(i)=c1;
            cerr2(i)=c2;
        end
    end
end

figure(); 
plot(freqs,coh);

if onlyDoDiags==0
    frs={'1';'2';'4';'6';'8';'10';'12';'14';'16';'18';'20';'30';'40';'50';'60'};
    figure();
    imagesc(coh_matrix);
    set(gca,'XTick',[1:size(coh_matrix,1)]);
    set(gca,'XTickLabel',frs);
    set(gca,'YTick',[1:size(coh_matrix,2)]);
    set(gca,'YTickLabel',frs);
    xlabel('Response Freqs (Hz)');
    ylabel('Stim Freqs (Hz)');
    title('Coherence Matrix');
    
    offDiags=zeros(1,size(coh_matrix,2));
    for i=1:size(coh_matrix,2)
        currRow=i;
        offDiags(i)=mean(coh_matrix([1:currRow-1 currRow+1:end],i));
    end
    offDiagMatrix=repmat(offDiags,size(coh_matrix,1),1);
    coh_matrix_minusOffDiags=coh_matrix-offDiagMatrix;
    figure();
    imagesc(coh_matrix_minusOffDiags);
    set(gca,'XTick',[1:size(coh_matrix,1)]);
    set(gca,'XTickLabel',frs);
    set(gca,'YTick',[1:size(coh_matrix,2)]);
    set(gca,'YTickLabel',frs);
    xlabel('Response Freqs (Hz)');
    ylabel('Stim Freqs (Hz)');
    title('Coherence Matrix Minus Off-Diagonals');
    
    diags=zeros(1,size(coh_matrix,2));
    for i=1:length(diags)
        diags(i)=coh_matrix_minusOffDiags(i,i);
    end
    figure();
    plot(freqs,diags);
    xlabel('Stim Freqs (Hz)');
    ylabel('Response');
end
