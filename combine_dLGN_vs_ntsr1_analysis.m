function combine_dLGN_vs_ntsr1_analysis(datadir)

placeStimOnsetAt=4; % seconds from trial onset

if iscell(datadir)
    
    for i=1:length(datadir)
        d=datadir{i};
        a=load([d '\' 'stimWindow.mat']);
        stimWindow=a.stimWindow;
        
        a=load([d '\dLGN\specs\' 'noTheta_lowF1_dLGN']);
        noTheta_lowF1_dLGN=a.takeLow;
        a=load([d '\dLGN\specs\' 'noTheta_highF1_dLGN']);
        noTheta_highF1_dLGN=a.takeHigh;
        a=load([d '\dLGN\specs\' 'theta_lowF1_dLGN']);
        theta_lowF1_dLGN=a.takeLow;
        a=load([d '\dLGN\specs\' 'theta_highF1_dLGN']);
        theta_highF1_dLGN=a.takeHigh;
        a=load([d '\dLGN\specs\' 'noTheta_noLED']);
        noTheta=a.noTheta;
        
        a=load([d '\V1 Ntsr1\specs\' 'noTheta_lowF1_Ntsr1']);
        noTheta_lowF1_Ntsr1=a.takeLow;
        a=load([d '\V1 Ntsr1\specs\' 'noTheta_highF1_Ntsr1']);
        noTheta_highF1_Ntsr1=a.takeHigh;
        a=load([d '\V1 Ntsr1\specs\' 'theta_lowF1_Ntsr1']);
        theta_lowF1_Ntsr1=a.takeLow;
        a=load([d '\V1 Ntsr1\specs\' 'theta_highF1_Ntsr1']);
        theta_highF1_Ntsr1=a.takeHigh;
        a=load([d '\dLGN\specs\' 'theta_noLED']);
        theta=a.theta;
        
        noTheta_lowF1_dLGN=alignToStim(noTheta_lowF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        noTheta_highF1_dLGN=alignToStim(noTheta_highF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        theta_lowF1_dLGN=alignToStim(theta_lowF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        theta_highF1_dLGN=alignToStim(theta_highF1_dLGN,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        
        noTheta_lowF1_Ntsr1=alignToStim(noTheta_lowF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        noTheta_highF1_Ntsr1=alignToStim(noTheta_highF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        theta_lowF1_Ntsr1=alignToStim(theta_lowF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        theta_highF1_Ntsr1=alignToStim(theta_highF1_Ntsr1,noTheta.allS.t,placeStimOnsetAt,stimWindow(1));
        
        if i==1
            all_noTheta_lowF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
            all_noTheta_highF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
            all_theta_lowF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
            all_theta_highF1_dLGN=zeros(size(noTheta_lowF1_dLGN));
            
            all_noTheta_lowF1_Ntsr1=zeros(size(noTheta_lowF1_dLGN));
            all_noTheta_highF1_Ntsr1=zeros(size(noTheta_lowF1_dLGN));
            all_theta_lowF1_Ntsr1=zeros(size(noTheta_lowF1_dLGN));
            all_theta_highF1_Ntsr1=zeros(size(noTheta_lowF1_dLGN));
        end
        
        
        all_noTheta_lowF1_dLGN=addAll(all_noTheta_lowF1_dLGN,length(noTheta.allS.S)*noTheta_lowF1_dLGN);
        all_noTheta_highF1_dLGN=addAll(all_noTheta_highF1_dLGN,length(noTheta.allS.S)*noTheta_highF1_dLGN);
        all_theta_lowF1_dLGN=addAll(all_theta_lowF1_dLGN,length(noTheta.allS.S)*theta_lowF1_dLGN);
        all_theta_highF1_dLGN=addAll(all_theta_highF1_dLGN,length(noTheta.allS.S)*theta_highF1_dLGN);
        
        all_noTheta_lowF1_Ntsr1=addAll(all_noTheta_lowF1_Ntsr1,length(noTheta.allS.S)*noTheta_lowF1_Ntsr1);
        all_noTheta_highF1_Ntsr1=addAll(all_noTheta_highF1_Ntsr1,length(noTheta.allS.S)*noTheta_highF1_Ntsr1);
        all_theta_lowF1_Ntsr1=addAll(all_theta_lowF1_Ntsr1,length(noTheta.allS.S)*theta_lowF1_Ntsr1);
        all_theta_highF1_Ntsr1=addAll(all_theta_highF1_Ntsr1,length(noTheta.allS.S)*theta_highF1_Ntsr1);
    end 
else
    error('expected datadir to be a cell array');
end

all_noTheta_lowF1_dLGN=all_noTheta_lowF1_dLGN(noTheta.allS.t<12,noTheta.allS.f<=50);
all_noTheta_highF1_dLGN=all_noTheta_highF1_dLGN(noTheta.allS.t<12,noTheta.allS.f<=50);
all_theta_lowF1_dLGN=all_theta_lowF1_dLGN(noTheta.allS.t<12,noTheta.allS.f<=50);
all_theta_highF1_dLGN=all_theta_highF1_dLGN(noTheta.allS.t<12,noTheta.allS.f<=50);

all_noTheta_lowF1_Ntsr1=all_noTheta_lowF1_Ntsr1(noTheta.allS.t<12,noTheta.allS.f<=50);
all_noTheta_highF1_Ntsr1=all_noTheta_highF1_Ntsr1(noTheta.allS.t<12,noTheta.allS.f<=50);
all_theta_lowF1_Ntsr1=all_theta_lowF1_Ntsr1(noTheta.allS.t<12,noTheta.allS.f<=50);
all_theta_highF1_Ntsr1=all_theta_highF1_Ntsr1(noTheta.allS.t<12,noTheta.allS.f<=50);

all_noTheta_lowF1_dLGN=whitenAndNormSpecgram(all_noTheta_lowF1_dLGN,false);
all_noTheta_highF1_dLGN=whitenAndNormSpecgram(all_noTheta_highF1_dLGN,false);
all_theta_lowF1_dLGN=whitenAndNormSpecgram(all_theta_lowF1_dLGN,false);
all_theta_highF1_dLGN=whitenAndNormSpecgram(all_theta_highF1_dLGN,false);

all_noTheta_lowF1_Ntsr1=whitenAndNormSpecgram(all_noTheta_lowF1_Ntsr1,true);
all_noTheta_highF1_Ntsr1=whitenAndNormSpecgram(all_noTheta_highF1_Ntsr1,true);
all_theta_lowF1_Ntsr1=whitenAndNormSpecgram(all_theta_lowF1_Ntsr1,true);
all_theta_highF1_Ntsr1=whitenAndNormSpecgram(all_theta_highF1_Ntsr1,true);


figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_lowF1_dLGN'); title('no theta dLGN low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_highF1_dLGN'); title('no theta dLGN high F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_lowF1_dLGN'); title('theta dLGN low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_highF1_dLGN'); title('theta dLGN high F1');

figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_lowF1_Ntsr1'); title('no theta Ntsr1 low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_noTheta_highF1_Ntsr1'); title('no theta Ntsr1 high F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_lowF1_Ntsr1'); title('theta Ntsr1 low F1');
figure(); imagesc(noTheta.allS.t(noTheta.allS.t<12),noTheta.allS.f(noTheta.allS.f<=50),all_theta_highF1_Ntsr1'); title('theta Ntsr1 high F1');

end

function temp=whitenAndNormSpecgram(temp,isCx)

if isCx==false
    temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
    temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
    temp(:,2)=1.07*temp(:,2);
    temp(:,4)=1.07*temp(:,4);
    temp(:,3)=1.0475*temp(:,3);
    % temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    % temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
else
    temp(:,1:3)=temp(:,1:3)./(3.3*0.17.*(0.7225/0.1121));
    temp(:,4)=temp(:,4)./(2.0*0.73.*(0.7225/0.5035));
    temp(:,2)=1.07*temp(:,2);
    temp(:,4)=1.07*temp(:,4);
%     temp(:,3)=1.0475*temp(:,3);
    temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
    temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
end

K=ones([2 2]); 
temp=conv2(temp,K,'same');
temp=interp2(temp,4);

temp=temp(1:end-13,1:end-13);

end

function out=addAll(data1,data2)

tmp=cat(3,data1,data2);
out=nansum(tmp,3);

end

function align_specgram=alignToStim(specgram,specgram_t,placeStimOnsetAt,currStimOnset)

if placeStimOnsetAt==currStimOnset
    align_specgram=specgram;
    return
end

timeDiff=abs(currStimOnset-placeStimOnsetAt);
inds_diff=floor(timeDiff/(specgram_t(2)-specgram_t(1)));

if currStimOnset>placeStimOnsetAt
    align_specgram=[specgram(inds_diff:end,:); nan(inds_diff-1,size(specgram,2))];
elseif currStimOnset<placeStimOnsetAt
    align_specgram=[nan(inds_diff,size(specgram,2)); specgram(1:end-inds_diff,:)];
end

end
       