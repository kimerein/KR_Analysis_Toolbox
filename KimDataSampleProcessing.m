% clear
% [~,name] = system('hostname');
% switch name(1:end-1)
%     case 'zpow'
        computerFolder = 'W:\Simultaneous dLGN and V1 for Narrowband Gamma';
%     otherwise  
%         display('Add case for computer name');
% end

processRecording = 315;

switch processRecording
    case 315
        % % % Recording 315
        % % % High Contrast
        folder = [computerFolder filesep 'Recording 315' filesep 'Multi-unit spikes and V1 LFP' filesep 'High contrast drifting grating'];
        spikeFile = [folder filesep 'dLGNandV1spikes'];
        LFP_file =  [folder filesep 'V1_LFP'];
        LFP_sampleRate_file = [folder filesep 'LFP_sampleRate'];
    case 312
        % % % Recording 312
        % % % High Contrast
        folder = [computerFolder filesep 'Recording 312' filesep 'Multi-unit spikes and V1 LFP' filesep 'High contrast drifting grating'];
        spikeFile = [folder filesep 'dLGNandV1spikes'];
        LFP_file =  [folder filesep 'V1_LFP'];
        LFP_sampleRate_file = [folder filesep 'LFP_sampleRate'];
    case 304
        % % % Recording 304
        % % % Low Contrast
        folder = [computerFolder filesep 'Recording 304' filesep 'Multi-unit spikes and V1 LFP' filesep 'Low contrast visual stimulus'];
        spikeFile = [folder filesep 'lowcontrast_dLGNandV1spikes'];
        LFP_file =  [folder filesep 'lowcontrast_V1_LFP'];
        LFP_sampleRate_file = [folder filesep 'LFP_sampleRate'];
end
display('Loading data...');
load(spikeFile);
load(LFP_file);
load(LFP_sampleRate_file);
clear spikeFile LFP_file LFP_sampleRate_file
% Reorder LGN spike times
if exist('dLGNspikes')
    clear LGN_spiketimes
    trialIDs = unique(dLGNspikes.trials);
    for itrial = 1:length(trialIDs)
        LGN_spiketimes{itrial} = dLGNspikes.spiketimes(dLGNspikes.trials == trialIDs(itrial));
    end
    LGN_spiketimes = cell2struct(LGN_spiketimes, 'times');
end
% Reorder V1 spike times
if exist('V1spikes')
    clear V1_spiketimes
    trialIDs = unique(V1spikes.trials);
    for itrial = 1:length(trialIDs)
        V1_spiketimes{itrial} = (V1spikes.spiketimes(V1spikes.trials == trialIDs(itrial)));
    end
    V1_spiketimes = cell2struct(V1_spiketimes, 'times');
end   

movingwin=[1 0.2];% [3 1];%[window_size window_shift] [3 1];
params.Fs = LFP_sampleRate;
params.tapers=[5 9];
params.fpass=[30 95];
params.trialave=0;

display('Getting coherence on V1 LFP and V1 Spiking...')
[vC,vphi,vSLfpSpikes,SLfp,vSSpikes,t,f]=cohgramcpt(LFPbySweep',V1_spiketimes,movingwin,params);
display('Getting coherence on V1 LFP and LGN Spiking...')
[lC,lphi,lSLfpSpikes,vSLfp,lSSpikes,t,f]=cohgramcpt(LFPbySweep',LGN_spiketimes,movingwin,params);
display('Done!')
%% Plot figure
figure(50)
subplot(2,2,1)
imagesc(t,f,log(nanmean(vSSpikes,3))');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
colorbar; axis xy; title('V1 Spikes')
ylabel('Freq (Hz)')
xlabel('Time')
subplot(2,2,2)
imagesc(t,f,log(nanmean(SLfp,3))');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
colorbar; axis xy; title('V1 LFP')
ylabel('Freq (Hz)')
xlabel('Time')
subplot(2,2,3)
imagesc(t,f,log(nanmean(lSSpikes,3))');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
colorbar; axis xy; title('LGN Spikes')
ylabel('Freq (Hz)')
xlabel('Time')
subplot(2,2,4)
imagesc(t,f,abs(nanmean(lSLfpSpikes,3))');
% imagesc(t,f,abs(nanmean(lC,3)));
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
colorbar; axis xy; title('V1 LFP - LGN Spikes')
ylabel('Freq (Hz)')
xlabel('Time')
% colormap