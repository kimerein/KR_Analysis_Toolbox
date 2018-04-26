function coherenceV1Thal_wrapper(V1_psths,thal_psths)

l_V1=length(V1_psths);
l_thal=length(thal_psths);
for i=1:l_V1
    for j=1:l_thal
        disp(i);
        disp(j);
        psths.V1_x=V1_psths(i).V1_x;
        psths.V1_y=V1_psths(i).V1_y;
        psths.thal_x=thal_psths(j).thal_x;
        psths.thal_y=thal_psths(j).thal_y;
        name=['thalV1coherence_' num2str(j) num2str(i)]; 
        coherenceV1PSTHwThalPSTH(name,[],[],psths);
    end
end