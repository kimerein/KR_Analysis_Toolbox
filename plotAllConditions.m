function [isSig_thetaFx,isSig_LEDnoTheta,isSig_LEDtheta]=plotAllConditions(noTheta_noLED,theta_noLED,noTheta_LED,theta_LED,bigFres,plotCell)

figure();
ha = tight_subplot(1,3,[.03 .05],[.05 .1],[.06 .03]);

allPlotLims=[];

axes(ha(1));
[plotLims,r1]=plotFreqResponseBootstrap(noTheta_noLED,bigFres,plotCell,'k');
allPlotLims=[allPlotLims plotLims];
[plotLims,r2]=plotFreqResponseBootstrap(theta_noLED,bigFres,plotCell,[0.7 0.5 0.5]);
allPlotLims=[allPlotLims plotLims];
isSig_thetaFx=nan(1,size(r1,2));
for i=1:size(r1,2)
    if all(isnan(r1(plotCell,i,:))) | all(isnan(r2(plotCell,i,:)))
        continue
    end
    isSig_thetaFx(i)=ranksum(reshape(r1(plotCell,i,:),1,length(r1(plotCell,i,:))),reshape(r2(plotCell,i,:),1,length(r2(plotCell,i,:))));
end

axes(ha(2));
[plotLims,r1]=plotFreqResponseBootstrap(noTheta_noLED,bigFres,plotCell,'k');
allPlotLims=[allPlotLims plotLims];
[plotLims,r2]=plotFreqResponseBootstrap(noTheta_LED,bigFres,plotCell,'b');
allPlotLims=[allPlotLims plotLims];
isSig_LEDnoTheta=nan(1,size(r1,2));
for i=1:size(r1,2)
    if all(isnan(r1(plotCell,i,:))) | all(isnan(r2(plotCell,i,:)))
        continue
    end
    isSig_LEDnoTheta(i)=ranksum(reshape(r1(plotCell,i,:),1,length(r1(plotCell,i,:))),reshape(r2(plotCell,i,:),1,length(r2(plotCell,i,:))));
end

axes(ha(3));
[plotLims,r1]=plotFreqResponseBootstrap(theta_noLED,bigFres,plotCell,[0.7 0.5 0.5]);
allPlotLims=[allPlotLims plotLims];
[plotLims,r2]=plotFreqResponseBootstrap(theta_LED,bigFres,plotCell,'c');
allPlotLims=[allPlotLims plotLims];
isSig_LEDtheta=nan(1,size(r1,2));
for i=1:size(r1,2)
    if all(isnan(r1(plotCell,i,:))) | all(isnan(r2(plotCell,i,:)))
        continue
    end
    isSig_LEDtheta(i)=ranksum(reshape(r1(plotCell,i,:),1,length(r1(plotCell,i,:))),reshape(r2(plotCell,i,:),1,length(r2(plotCell,i,:))));
end

axes(ha(1));
xlim([1 60]);
ylim([nanmin(allPlotLims) nanmax(allPlotLims)]);

axes(ha(2));
xlim([1 60]);
ylim([nanmin(allPlotLims) nanmax(allPlotLims)]);

axes(ha(3));
xlim([1 60]);
ylim([nanmin(allPlotLims) nanmax(allPlotLims)]);

set(gcf, 'Position', [100, 500, 700*0.8, 300*0.8]);

end

