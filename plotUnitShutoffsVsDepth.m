function lastNonzero=plotUnitShutoffsVsDepth(spikes,fileInd)

ledOnset=1.34;
dropTo=0.3679;
until=1;

a=unique(spikes.assigns);
lastNonzero=zeros(1,length(a));
for i=1:length(a)
    [xpoints,ypoints1,ypoints2]=scriptForComparingMUA(filtspikes(spikes,0,'assigns',a(i)),fileInd,[],[]);
    ypoints2=smooth(ypoints2,10);
    ypoints2=smooth(ypoints2,20);
%     ypoints2=smooth(ypoints2,3);
%     ypoints2=smooth(ypoints2,5);
    subx=xpoints(xpoints>=ledOnset-0.01 & xpoints<=ledOnset+until);
%     nonz=find(ypoints2(xpoints>=ledOnset-0.01 & xpoints<=ledOnset+until)>0);
    indAtLedOnset=find(xpoints<ledOnset,1,'last');
    peakVal=ypoints2(indAtLedOnset);
    baseVal=mean(ypoints2(xpoints>=0.7 & xpoints<0.9));
    threshVal=(peakVal-baseVal)*dropTo;
    nonz=find(ypoints2(xpoints>=ledOnset-0.01 & xpoints<=ledOnset+until)<threshVal,1,'first');
    if isempty(nonz)
        lastNonzero(i)=subx(end);
    else
        lastNonzero(i)=subx(nonz(end));
    end
    if i==1
        figure(); 
        plot(xpoints,ypoints2,'Color','r');
        hold on; 
        plot([xpoints(1) xpoints(end)],[threshVal threshVal]);
        scatter(lastNonzero(i),threshVal);
    end
end

lastNonzero=lastNonzero-ledOnset;