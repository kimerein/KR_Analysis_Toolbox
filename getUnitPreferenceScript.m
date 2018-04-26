tryS=1:12;

psth=psthGp1;
autocorr_out=autocorr_outGp1;

sResponses=zeros(1,length(tryS));
for z=1:length(sResponses)
    uses=tryS(z);
    scriptForPlottingAutocorr;
    close all
    sResponses(z)=outS;
end

figure(); 
plot(tryS,sResponses);