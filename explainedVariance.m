function explainedVariance(dLGNpsth,V1psth)

doTrialAve=0;
useStimInstead=0;

useStimWindow=500:9600; % in inds
dLGNpsth=dLGNpsth(:,useStimWindow);
V1psth=V1psth(:,useStimWindow);

bins=[1:9 10:10:90 100:10:1000];

explainedVariance=zeros(size(bins));
for i=1:length(bins)
    currBin=bins(i);
    if doTrialAve==1
        a=downSampAv(nanmean(dLGNpsth,1),currBin)';
        b=downSampAv(nanmean(V1psth,1),currBin)';
        if useStimInstead==1
            x=linspace(0,3,floor(length(a)/3));
            c=chirp(x,1,3,30,'logarithmic')+0.5;
            b=[c fliplr(c) c];
            b=[b ones(1,length(a)-length(b)).*c(end)]';
        end
    else
        a_temp=downSampMatrix(dLGNpsth,currBin);
        b_temp=downSampMatrix(V1psth,currBin);
        a=a_temp(1:end)';
        b=b_temp(1:end)';
    end
    r=corrcoef([a b]);
%     r2=corrcoef([downSampAv(nanmean(dLGNpsth(floor(currBin/4)+1:end,:),1),currBin)' downSampAv(nanmean(V1psth(floor(currBin/4)+1:end,:),1),currBin)']);
%     r3=corrcoef([downSampAv(nanmean(dLGNpsth(floor(2*currBin/4)+1:end,:),1),currBin)' downSampAv(nanmean(V1psth(floor(2*currBin/4)+1:end,:),1),currBin)']);
%     r4=corrcoef([downSampAv(nanmean(dLGNpsth(floor(3*currBin/4)+1:end,:),1),currBin)' downSampAv(nanmean(V1psth(floor(3*currBin/4)+1:end,:),1),currBin)']);
%     allr=nanmean([r(1,2) r2(1,2) r3(1,2) r4(1,2)]);
    explainedVariance(i)=r(1,2)^2;
    if currBin==900
        figure(); 
        scatter(a+100*rand(size(a)),b+100*rand(size(b)));
%         scatter(a,b);
    end
end

figure(); 
plot(bins,explainedVariance);
xlabel('Time Bin Size (ms)');
ylabel('Explained Variance');

figure(); 
semilogx(1./(bins./1000),explainedVariance);
xlabel('Frequency (Hz)');
ylabel('Explained Variance');