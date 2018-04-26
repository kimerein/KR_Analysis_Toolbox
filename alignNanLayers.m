function [normedLayers,meanAcrossRows]=alignNanLayers(allNanLayers,useRows,freqs)

meanAcrossRows=nanmean(allNanLayers(useRows,:,:),1);
maxForEachMouse=nanmax(meanAcrossRows,[],2);
normedLayers=meanAcrossRows;
for j=1:size(meanAcrossRows,3)
    for i=1:length(freqs)
        normedLayers(1,i,j)=normedLayers(1,i,j)/maxForEachMouse(1,1,j);
    end
end

figure(); 
semilogx(freqs,nanmean(normedLayers,3));
hold on; 
semilogx(freqs,nanmean(normedLayers,3)+nanstd(normedLayers,[],3)./sqrt(size(normedLayers,3)));
semilogx(freqs,nanmean(normedLayers,3)-nanstd(normedLayers,[],3)./sqrt(size(normedLayers,3)));