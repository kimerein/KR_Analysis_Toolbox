function trialByTrialFreqCorr(dd)

% for each pair of frequency bands, calculate the correlation coefficient during stim vs. spont
freqbands=[0 2; 2 4; 4 6; 6 8; 8 10; 10 12; 12 14; 14 16; 16 18; 18 20; 20 22; 22 24; 24 26; 26 28; 28 30; 30 32; ...
           32 34; 34 36; 36 38; 38 40; 40 42; 42 44; 44 46; 46 48; 48 49.999];

ls=dir(dd);
rhos_stim_noThetaNoLED=nan(size(freqbands,1),size(freqbands,1),1);
rhos_spont_noThetaNoLED=nan(size(freqbands,1),size(freqbands,1),1);
whichcell_noThetaNoLED=1;
rhos_stim_noThetaLED=nan(size(freqbands,1),size(freqbands,1),1);
rhos_spont_noThetaLED=nan(size(freqbands,1),size(freqbands,1),1);
whichcell_noThetaLED=1;
rhos_stim_thetaNoLED=nan(size(freqbands,1),size(freqbands,1),1);
rhos_spont_thetaNoLED=nan(size(freqbands,1),size(freqbands,1),1);
whichcell_thetaNoLED=1;
rhos_stim_thetaLED=nan(size(freqbands,1),size(freqbands,1),1);
rhos_spont_thetaLED=nan(size(freqbands,1),size(freqbands,1),1);
whichcell_thetaLED=1;
for i=1:length(ls)
    na=ls(i).name;
    if ~isempty(regexp(na(1),'M'))
        % data directory
        % load stim window
        if ~exist([ls(i).folder '\' ls(i).name '\stimWindow.mat'])
            continue
        end
        if ~exist([ls(i).folder '\' ls(i).name '\dLGN'])
            continue
        end
        if ~exist([ls(i).folder '\' ls(i).name '\V1 Ntsr1'])
            continue
        end
        a=load([ls(i).folder '\' ls(i).name '\stimWindow.mat']);
        stimWindow=a.stimWindow;
        % load t and f for dLGN data
        a=load([ls(i).folder '\' ls(i).name '\dLGN\specs\noTheta_LED\allS.mat']);
        dLGN_t=a.temp.t;
        dLGN_f=a.temp.f;
        % load t and f for V1 Nstr1 data
        a=load([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\noTheta_LED\allS.mat']);
        V1_t=a.temp.t;
        V1_f=a.temp.f;
        % get freqinds
        dLGN_freqinds=nan(size(freqbands));
        V1_freqinds=nan(size(freqbands));
        for j=1:length(freqbands)
            currband=freqbands(j,:);
            dLGN_freqinds(j,:)=[find(dLGN_f>=currband(1),1,'first') find(dLGN_f>=currband(2),1,'first')];
            V1_freqinds(j,:)=[find(V1_f>=currband(1),1,'first') find(V1_f>=currband(2),1,'first')];
        end
        % get inds for spont v. stim
        dLGN_stim_inds=find(dLGN_t>stimWindow(1) & dLGN_t<=stimWindow(2));
        V1_stim_inds=find(V1_t>stimWindow(1) & V1_t<=stimWindow(2));
        dLGN_2nd_spont_inds=find(dLGN_t>stimWindow(2)+1);
        V1_2nd_spont_inds=find(V1_t>stimWindow(2)+1);
        % take smaller size
        if dLGN_stim_inds(end)<V1_stim_inds(end)
            consensus_stim_inds=dLGN_stim_inds;
        elseif dLGN_stim_inds(end)>=V1_stim_inds(end)
            consensus_stim_inds=V1_stim_inds;
        end
        if dLGN_2nd_spont_inds(1)<V1_2nd_spont_inds(1)
            consensus_2nd_spont_inds=dLGN_2nd_spont_inds;
        else
            consensus_2nd_spont_inds=V1_2nd_spont_inds;
        end
        %consensus_spont_inds=[1:consensus_stim_inds(1)-1 consensus_2nd_spont_inds];
        consensus_spont_inds=[1:consensus_stim_inds(1)-1]; % take only pre-stim
        
        % No theta no LED
        disp('No theta no LED');
        for j=1:(length(dir([ls(i).folder '\' ls(i).name '\dLGN\specs\noTheta_noLED_allS_S']))-2)
            % load dLGN cell by cell
            disp(j);
            a=load([ls(i).folder '\' ls(i).name '\dLGN\specs\noTheta_noLED_allS_S\' 'cell' num2str(j) '.mat']);
            dLGN_cell=a.temp;
            for k=1:(length(dir([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\noTheta_noLED_allS_S']))-2)
                % load V1 Ntsr1 cell by cell
                a=load([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\noTheta_noLED_allS_S\' 'cell' num2str(k) '.mat']);
                V1_cell=a.temp;
                % for each pair of frequency bands, during spont v. stim, calculate
                % rho, corr coeff
                if isempty(dLGN_cell) || isempty(V1_cell)
                    rhos_stim_noThetaNoLED(size(freqbands,1),size(freqbands,1),whichcell_noThetaNoLED)=nan;
                    rhos_spont_noThetaNoLED(size(freqbands,1),size(freqbands,1),whichcell_noThetaNoLED)=nan;
                    whichcell_noThetaNoLED=whichcell_noThetaNoLED+1;
                    continue
                end
                for f1=1:size(freqbands,1)
                    for f2=1:size(freqbands,1)
                        % stim
                        % first dim is time
                        % second dim is freq
                        % third dim is trial
                        currband1=nanmean(dLGN_cell(consensus_stim_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
                        currband1=currband1(1:end);
                        currband2=nanmean(V1_cell(consensus_stim_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
                        currband2=currband2(1:end);
                        rho=corr(currband1',currband2');
                        rhos_stim_noThetaNoLED(f1,f2,whichcell_noThetaNoLED)=rho;
                        currband1=nanmean(dLGN_cell(consensus_spont_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
                        currband1=currband1(1:end);
                        currband2=nanmean(V1_cell(consensus_spont_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
                        currband2=currband2(1:end);
                        rho=corr(currband1',currband2');
                        rhos_spont_noThetaNoLED(f1,f2,whichcell_noThetaNoLED)=rho;
                    end
                end
                whichcell_noThetaNoLED=whichcell_noThetaNoLED+1;
            end
        end
        
        % No theta LED
%         disp('No theta LED');
%         for j=1:(length(dir([ls(i).folder '\' ls(i).name '\dLGN\specs\noTheta_LED_allS_S']))-2)
%             % load dLGN cell by cell
%             disp(j);
%             a=load([ls(i).folder '\' ls(i).name '\dLGN\specs\noTheta_LED_allS_S\' 'cell' num2str(j) '.mat']);
%             dLGN_cell=a.temp;
%             for k=1:(length(dir([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\noTheta_LED_allS_S']))-2)
%                 % load V1 Ntsr1 cell by cell
%                 a=load([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\noTheta_LED_allS_S\' 'cell' num2str(k) '.mat']);
%                 V1_cell=a.temp;
%                 % for each pair of frequency bands, during spont v. stim, calculate
%                 % rho, corr coeff
%                 if isempty(dLGN_cell) || isempty(V1_cell)
%                     rhos_stim_noThetaLED(size(freqbands,1),size(freqbands,1),whichcell_noThetaLED)=nan;
%                     rhos_spont_noThetaLED(size(freqbands,1),size(freqbands,1),whichcell_noThetaLED)=nan;
%                     whichcell_noThetaLED=whichcell_noThetaLED+1;
%                     continue
%                 end
%                 for f1=1:size(freqbands,1)
%                     for f2=1:size(freqbands,1)
%                         % stim
%                         % first dim is time
%                         % second dim is freq
%                         % third dim is trial
%                         currband1=nanmean(dLGN_cell(consensus_stim_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
%                         currband1=currband1(1:end);
%                         currband2=nanmean(V1_cell(consensus_stim_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
%                         currband2=currband2(1:end);
%                         rho=corr(currband1',currband2');
%                         rhos_stim_noThetaLED(f1,f2,whichcell_noThetaLED)=rho;
%                         currband1=nanmean(dLGN_cell(consensus_spont_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
%                         currband1=currband1(1:end);
%                         currband2=nanmean(V1_cell(consensus_spont_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
%                         currband2=currband2(1:end);
%                         rho=corr(currband1',currband2');
%                         rhos_spont_noThetaLED(f1,f2,whichcell_noThetaLED)=rho;
%                     end
%                 end
%                 whichcell_noThetaLED=whichcell_noThetaLED+1;
%             end
%         end
        
        % Theta no LED
        disp('Theta no LED');
        for j=1:(length(dir([ls(i).folder '\' ls(i).name '\dLGN\specs\theta_noLED_allS_S']))-2)
            % load dLGN cell by cell
            disp(j);
            a=load([ls(i).folder '\' ls(i).name '\dLGN\specs\theta_noLED_allS_S\' 'cell' num2str(j) '.mat']);
            dLGN_cell=a.temp;
            for k=1:(length(dir([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\theta_noLED_allS_S']))-2)
                % load V1 Ntsr1 cell by cell
                a=load([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\theta_noLED_allS_S\' 'cell' num2str(k) '.mat']);
                V1_cell=a.temp;
                % for each pair of frequency bands, during spont v. stim, calculate
                % rho, corr coeff
                if isempty(dLGN_cell) || isempty(V1_cell)
                    rhos_stim_thetaNoLED(size(freqbands,1),size(freqbands,1),whichcell_thetaNoLED)=nan;
                    rhos_spont_thetaNoLED(size(freqbands,1),size(freqbands,1),whichcell_thetaNoLED)=nan;
                    whichcell_thetaNoLED=whichcell_thetaNoLED+1;
                    continue
                end
                for f1=1:size(freqbands,1)
                    for f2=1:size(freqbands,1)
                        % stim
                        % first dim is time
                        % second dim is freq
                        % third dim is trial
                        currband1=nanmean(dLGN_cell(consensus_stim_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
                        currband1=currband1(1:end);
                        currband2=nanmean(V1_cell(consensus_stim_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
                        currband2=currband2(1:end);
                        rho=corr(currband1',currband2');
                        rhos_stim_thetaNoLED(f1,f2,whichcell_thetaNoLED)=rho;
                        currband1=nanmean(dLGN_cell(consensus_spont_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
                        currband1=currband1(1:end);
                        currband2=nanmean(V1_cell(consensus_spont_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
                        currband2=currband2(1:end);
                        rho=corr(currband1',currband2');
                        rhos_spont_thetaNoLED(f1,f2,whichcell_thetaNoLED)=rho;
                    end
                end
                whichcell_thetaNoLED=whichcell_thetaNoLED+1;
            end
        end
        
        % Theta LED
%         disp('Theta LED');
%         for j=1:(length(dir([ls(i).folder '\' ls(i).name '\dLGN\specs\theta_LED_allS_S']))-2)
%             % load dLGN cell by cell
%             disp(j);
%             a=load([ls(i).folder '\' ls(i).name '\dLGN\specs\theta_LED_allS_S\' 'cell' num2str(j) '.mat']);
%             dLGN_cell=a.temp;
%             for k=1:(length(dir([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\theta_LED_allS_S']))-2)
%                 % load V1 Ntsr1 cell by cell
%                 a=load([ls(i).folder '\' ls(i).name '\V1 Ntsr1\specs\theta_LED_allS_S\' 'cell' num2str(k) '.mat']);
%                 V1_cell=a.temp;
%                 % for each pair of frequency bands, during spont v. stim, calculate
%                 % rho, corr coeff
%                 if isempty(dLGN_cell) || isempty(V1_cell)
%                     rhos_stim_thetaLED(size(freqbands,1),size(freqbands,1),whichcell_thetaLED)=nan;
%                     rhos_spont_thetaLED(size(freqbands,1),size(freqbands,1),whichcell_thetaLED)=nan;
%                     whichcell_thetaLED=whichcell_thetaLED+1;
%                     continue
%                 end
%                 for f1=1:size(freqbands,1)
%                     for f2=1:size(freqbands,1)
%                         % stim
%                         % first dim is time
%                         % second dim is freq
%                         % third dim is trial
%                         currband1=nanmean(dLGN_cell(consensus_stim_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
%                         currband1=currband1(1:end);
%                         currband2=nanmean(V1_cell(consensus_stim_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
%                         currband2=currband2(1:end);
%                         rho=corr(currband1',currband2');
%                         rhos_stim_thetaLED(f1,f2,whichcell_thetaLED)=rho;
%                         currband1=nanmean(dLGN_cell(consensus_spont_inds,dLGN_freqinds(f1,1):dLGN_freqinds(f1,2),:),2);
%                         currband1=currband1(1:end);
%                         currband2=nanmean(V1_cell(consensus_spont_inds,V1_freqinds(f2,1):V1_freqinds(f2,2),:),2);
%                         currband2=currband2(1:end);
%                         rho=corr(currband1',currband2');
%                         rhos_spont_thetaLED(f1,f2,whichcell_thetaLED)=rho;
%                     end
%                 end
%                 whichcell_thetaLED=whichcell_thetaLED+1;
%             end
%         end
    end
end

plotResults(rhos_stim_noThetaNoLED,rhos_spont_noThetaNoLED,freqbands,'No theta no LED');
% plotResults(rhos_stim_noThetaLED,rhos_spont_noThetaLED,freqbands,'No theta LED');
plotResults(rhos_stim_thetaNoLED,rhos_spont_thetaNoLED,freqbands,'Theta no LED');
% plotResults(rhos_stim_thetaLED,rhos_spont_thetaLED,freqbands,'Theta LED');

end
            
function plotResults(rhos_stim_noThetaNoLED,rhos_spont_noThetaNoLED,freqbands,tit)

% adjust so that mean between pairs is zero correlation
takefornow=rhos_stim_noThetaNoLED;
tosubtract=nan(1,size(takefornow,3));
for dim1=1:size(takefornow,3)
    tosubtract(dim1)=nanmean(nanmean(takefornow(:,:,dim1),1),2);
end
for dim1=1:size(takefornow,3)
    takefornow(:,:,dim1)=takefornow(:,:,dim1)-tosubtract(dim1);
end
figure(); 
imagesc(nanmean(freqbands,2),nanmean(freqbands,2),reshape(nanmean(rhos_stim_noThetaNoLED(:,:,:),3),25,25)); 
ylabel('dLGN freq'); xlabel('V1 freq');
title([tit ' stim average']);
figure(); 
imagesc(nanmean(freqbands,2),nanmean(freqbands,2),reshape(nanmean(rhos_spont_noThetaNoLED(:,:,:),3),25,25)); 
ylabel('dLGN freq'); xlabel('V1 freq');
title([tit ' spont average']);
% Plot correlation after removing V1 freq power
takefornow=reshape(nanmean(rhos_stim_noThetaNoLED(:,:,:),3),25,25);
takefornow=takefornow./repmat(nanmean(takefornow,1),25,1);
figure();
imagesc(nanmean(freqbands,2),nanmean(freqbands,2),takefornow); 
ylabel('dLGN freq'); xlabel('V1 freq');
title([tit ' stim average -- norm by total V1 correlation per band']);
% Plot stim minus spont
takefornow=reshape(nanmean(rhos_stim_noThetaNoLED-rhos_spont_noThetaNoLED,3),25,25);
figure();
imagesc(nanmean(freqbands,2),nanmean(freqbands,2),takefornow); 
ylabel('dLGN freq'); xlabel('V1 freq');
title([tit ' stim minus spont']);

end
                
                