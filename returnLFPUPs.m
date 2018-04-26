function UPstates_LFP=returnLFPUPs(xpoints,powerRatio,thresh)

UPstates_LFP=cell(length(powerRatio),1);
for i=1:length(powerRatio)
    UPstates_LFP{i}=returnUPs(xpoints,powerRatio{i},thresh);
end