function getNewUnitsFx(bestSpikes,spikes,newWindow,params,saveDirName)

if isempty(newWindow)
%     newWindow=[0 0.999];
    newWindow=[1.3 3];
%     newWindow=[2 2.999];
%     newWindow=[1.34 1.8];
end
[~,led1Spikes,led2Spikes,ns1,ns2,ass]=getUnitByUnitEffects(bestSpikes,spikes,unique(bestSpikes.assigns),params.useLEDcond{1},params.useLEDcond{2},params.useStimcond,params.useStimcond,newWindow);
unitsFx.led1Spikes=led1Spikes;
unitsFx.led2Spikes=led2Spikes;
unitsFx.n1=ns1;
unitsFx.n2=ns2;
unitsFx.assigns=ass;
save([saveDirName '\baseUnitsFx.mat'],'unitsFx');
% save([saveDirName '\unitsFx.mat'],'unitsFx');
    