function plotSingleTrialsFromPSTH(psth,useTrials)

ds=3;

tri=psth.unitTrials{1};
l=psth.unitLED{1};
t=psth.t;

summedP=zeros(size(psth.psths{1}));
for i=1:length(psth.psths)
    p=psth.psths{i};
    summedP=summedP+p;
end

figure();
runningY=0;
hold on;
if ~isempty(ds)
    t=downSampAv(t,ds);
end
for i=1:length(useTrials)
    curry=summedP(tri==useTrials(i),:);
    if ~isempty(ds)
        curry=downSampAv(curry,ds);
    end
    maxy=max(curry);  
    plot(t,curry+runningY,'Color','k');
    runningY=runningY+maxy;
end
    
figure(); 
plot(t,downSampAv(nanmean(summedP(ismember(tri,useTrials),:),1),ds));

        