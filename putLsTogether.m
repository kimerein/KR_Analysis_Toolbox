function newL=putLsTogether(l1,l2,l1_block,l2_block)

for j=1:length(l1)
    for i=1:length(l1{1})
        curr1=l1{j};
        curr2=l2{j};
        oknow{i}.FRs_noLED{1}=curr1{i}.FRs_LED{l1_block};
        oknow{i}.FRs_LED{1}=curr2{i}.FRs_LED{l2_block};
    end
    newL{j}=oknow;
end