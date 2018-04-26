function identifyFSunits(useDir,listing)

f=figure(); 
movegui(f,'south');
for i=1:length(listing)
    currFolder=listing(i).name;
%     if ~exist([useDir '\' currFolder '\bestUnitsInfo.mat'],'file')
    if ~exist([useDir '\' currFolder '\criterionUnitsInfo.mat'],'file')
        continue
    end
%     a=load([useDir '\' currFolder '\bestUnitsInfo.mat']);
    a=load([useDir '\' currFolder '\criterionUnitsInfo.mat']);
    newassignsinfo=a.newassignsinfo;
    newassignsinfo.isFS=nan(1,length(newassignsinfo.waveformWidths));
    for j=1:size(newassignsinfo.waveforms,1)
        plot(newassignsinfo.waveforms(j,:));
        isFS=questdlg('Is unit fast-spiking (FS)?');
        if strcmp(isFS,'No')
            newassignsinfo.isFS(j)=0;
        elseif strcmp(isFS,'Yes')
            newassignsinfo.isFS(j)=1;
        elseif strcmp(isFS,'Cancel')
            return
        end
    end
%     save([useDir '\' currFolder '\bestUnitsInfo.mat'],'newassignsinfo');
    save([useDir '\' currFolder '\criterionUnitsInfo.mat'],'newassignsinfo');
end