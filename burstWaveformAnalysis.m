function burstWaveformAnalysis(datadir)

if ~iscell(datadir)
    temp{1}=datadir;
    datadir=temp;
end

waveforms_x=[];
waveforms_y=[];
normwaveforms=[];

for i=1:length(datadir)
    d=datadir{i};
    a=load([d '\' 'dLGNspikes_best']);
    dLGNspikes=a.dLGNspikes;
    
    useAssigns=unique(dLGNspikes.assigns);
    [ma]=max(dLGNspikes.waveforms,[],2);
    [~,mi]=max(ma,[],3);
    for j=1:length(useAssigns)
        subspikes=filtspikes(dLGNspikes,0,'assigns',useAssigns(j));
        
        bestEvCh=mode(mi(ismember(dLGNspikes.assigns,useAssigns(j))));
        wvfmAmps=reshape(nanmin(subspikes.waveforms(:,:,bestEvCh),[],2),size(subspikes.waveforms,1),1);
        [n,xout]=hist(wvfmAmps,500);
        waveforms_x{
        plot(xout,n,'Color','k');