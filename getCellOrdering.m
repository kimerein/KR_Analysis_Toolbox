function [responses,baselines,evNoSubBase,responses_tbyt,baselines_tbyt,justt_resp,justt_base]=getCellOrdering(psth)

s=1;
stimWindow=[4 6.5];
% stimWindow=[1.5 5];
% baseWindow=[0 1.5];
baseWindow=[0 4];
baseSubtract=1;
getTrialByTrial=0;
usel=[5.05];

% takeTrials=[13:72 141:212 277:316]-12;
% if ~isempty(takeTrials)
%     newpsth.t=psth.t;
%     for i=1:length(psth.psths)
%         p=psth.psths{i};
%         newpsth.psths{i}=p(takeTrials,:);
%         p=psth.unitTrials{i};
%         newpsth.unitTrials{i}=p(takeTrials);
%         p=psth.unitStimcond{i};
%         newpsth.unitStimcond{i}=p(takeTrials);
%         p=psth.unitLED{i};
%         newpsth.unitLED{i}=p(takeTrials);
%     end
%     psth=newpsth;
% end



responses=zeros(length(psth),length(s));
baselines=zeros(length(psth),length(s));
for i=1:length(psth.psths)
    p=psth.psths{i};
    uTrials=psth.unitTrials{i};
    uStimcond=psth.unitStimcond{i};
    uLED=psth.unitLED{i};
    for j=1:length(s)
        responses(i,j)=nanmean(nanmean(p(uStimcond==s(j) & uLED==usel,psth.t>=stimWindow(1) & psth.t<=stimWindow(2)),1),2);
        baselines(i,j)=nanmean(nanmean(p(uStimcond==s(j) & uLED==usel,psth.t>=baseWindow(1) & psth.t<=baseWindow(2)),1),2);
        if getTrialByTrial==1
            responses_tbyt{i,j}=nanmean(p(uStimcond==s(j) & uLED==usel,psth.t>=stimWindow(1) & psth.t<=stimWindow(2)),2);
            baselines_tbyt{i,j}=nanmean(p(uStimcond==s(j) & uLED==usel,psth.t>=baseWindow(1) & psth.t<=baseWindow(2)),2);
        end
    end
end

if getTrialByTrial==1
    for i=1:length(psth.psths)
        te=psth.psths{i};
        justt_resp{i}=nanmean(te(:,psth.t>=stimWindow(1) & psth.t<=stimWindow(2)),2);
        justt_base{i}=nanmean(te(:,psth.t>=baseWindow(1) & psth.t<=baseWindow(2)),2);
    end
end

if baseSubtract==1
    evNoSubBase=responses;
    responses=responses-baselines;
end


    












% function getCellOrdering(psth,useUnit,refFreqBand,useWindow,useWindowForAmp,ampOrdering,F1phases_input,usel)
% 
% % takeN=3;
% % stimconds=1:12;
% % usel=[0.05];
% takeN=5;
% stimconds=1:12;
% % takeN=12;
% % stimconds=1:16;
% % usel=[0];
% showFigs=0;
% 
% psths=psth.psths{useUnit};
% trialLED=psth.unitLED{useUnit};
% trialStimcond=psth.unitStimcond{useUnit};
% trials=psth.unitTrials{useUnit};
% 
% params.Fs=1/(psth.t(2)-psth.t(1));
% params.tapers=[2 5];
% params.trialave=1;
% % params.fpass=[2.5 30];
% params.fpass=[1 30];
% 
% F1amps=nan(length(stimconds),1);
% F1phases=nan(length(stimconds),1);
% F1coherence=nan(length(stimconds),1);
% subt=psth.t(psth.t>=useWindowForAmp(1) & psth.t<=useWindowForAmp(2));
% refx=subt;
% refy=sin(2*pi*nanmean(refFreqBand).*refx);
% for i=1:length(stimconds)
%     currs=stimconds(i);
%     useTrials=ismember(trialLED,usel) & ismember(trialStimcond,currs);
%     temp=psths(useTrials,psth.t>=useWindowForAmp(1) & psth.t<=useWindowForAmp(2));
%     % Get F1 amp
%     [S,f]=mtspectrumpb(temp',params);
% %     [S,f]=mtspectrumpb((temp-repmat(nanmean(temp,2),1,size(temp,2)))',params);
% %     F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)))-nanmean([nanmean(S(f<refFreqBand(1))) nanmean(S(f>refFreqBand(2)))]);
% %     F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)))./nanmean(S);
% %     DCchange=nanmean(nanmean(temp,1),2);
% %     F1amps(i)=(sqrt(nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2))))/25.05)-DCchange;
%     F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)));
%     % Get F1 phase
%     [C,phi,~,~,~,f]=coherencycpb(repmat(refy',1,size(temp,1)),temp',params);
%     F1coherence(i)=nanmean(C(f>=refFreqBand(1) & f<=refFreqBand(2)));
%     F1phases(i)=nanmean(phi(f>=refFreqBand(1) & f<=refFreqBand(2)));
% end
% 
% [~,ampOrder]=sort(F1amps);
% if ~isempty(ampOrdering)
%     ampOrder=ampOrdering;
% end
% if ~isempty(F1phases_input)
%     F1phases=F1phases_input;
% end
% % takeTop=find(F1phases>=-pi/2 & F1phases<0);
% takeTop=ampOrder(end-takeN+1:end);