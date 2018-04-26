function [finterp,allSfirst,allSsecond,ffr]=getAutocorrSpectraScript(psthFirst,autocorr_outFirst,psthSecond,autocorr_outSecond)

% Get unit preferences
tryS=1:12;

psth=psthSecond;
autocorr_out=autocorr_outSecond;
finterp=3:0.1:20-0.1;
lagsWindow=[0 1500];
usel=[0];
u=1:length(psth.psths);
ffr=zeros(length(u),length(tryS));
for i=1:length(u)
    for z=1:length(tryS)
        uses=tryS(z);
        [f,S]=scriptForPlottingAutocorr(psth,autocorr_out,usel,uses,i,lagsWindow);
        close all
        currS=interp1(f,S,finterp);
        autocorr_ffresponse(z,:)=currS;
    end
    ffr(i,:)=nanmean(autocorr_ffresponse(:,finterp>=3 & finterp<=20),2);
end
ffr=ffr(:,[12 1:11]);

tryS=1:12;

psth=psthFirst;
u=1:length(psth.psths);
allSfirst=cell(1,length(u));
allSsecond=cell(1,length(u));
for i=1:length(u)
    psth=psthFirst;
    autocorr_out=autocorr_outFirst;
%     psth=psthSecond;
%     autocorr_out=autocorr_outSecond;
    finterp=3:0.1:20-0.1;
    clear autocorr_spectraFirst
    lagsWindow=[0 3000];
    usel=[0];
    for z=1:length(tryS)
        uses=tryS(z);
%         [f,S]=scriptForPlottingAutocorr(psth,autocorr_out,usel,1:12,i,lagsWindow);
        [f,S]=scriptForPlottingAutocorr(psth,autocorr_out,usel,uses,i,lagsWindow);
        close all  
        currS=interp1(f,S,finterp);
        autocorr_spectraFirst(z,:)=currS;
    end
    allSfirst{i}=autocorr_spectraFirst;
    
    psth=psthSecond;
    autocorr_out=autocorr_outSecond;
%     psth=psthFirst;
%     autocorr_out=autocorr_outFirst;
    clear autocorr_spectraSecond
    lagsWindow=[750 1500];
    usel=[0];
    for z=1:length(tryS)
        uses=tryS(z);
        [f,S]=scriptForPlottingAutocorr(psth,autocorr_out,usel,uses,i,lagsWindow);
        close all
        currS=interp1(f,S,finterp);
        autocorr_spectraSecond(z,:)=currS;
    end
    allSsecond{i}=autocorr_spectraSecond;
    % clear allSpectra
    % for z=1:length(tryS)
    %     if z+1>length(tryS)
    %         allSpectra(z,:)=nanmean([autocorr_spectraFirst(1,:); autocorr_spectraSecond(12,:)],1);
    %     else
    %         allSpectra(z,:)=nanmean([autocorr_spectraFirst(z+1,:); autocorr_spectraSecond(z,:)],1);
    %     end
    % end
end
    