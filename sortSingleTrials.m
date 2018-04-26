function trialClassifications=sortSingleTrials(bandPassedLFPbySweep,ledForSweeps,LFP_Fs,baselineWindow,precedingWindow,precedingGammaThresh,duringWindow,duringGammaThresh)

% precedingWindow is the time window preceding the LED pulse to examine for
% high vs. low gamma
% duringWindow is the time window during and just after the LED pulse to examine for
% high vs. low gamma
%
% Function returns trialClassifications, which is vector with a classification for
% each trial:
% Classification    LED      Gamma preceding        Gamma during
% 1                 off      low                    high
% 2                 off      low                    low
% 3                 off      high                   high
% 4                 off      high                   low
% 5                 on       low                    high
% 6                 on       low                    low
% 7                 on       high                   high
% 8                 on       high                   low

% precedingGammaThresh=1; % This number times the stdev of baseline gamma power proxy
%                         % plus the average baseline gamma power proxy 
%                         % will decide low vs. high gamma preceding LED
%                         % pulse
% duringGammaThresh=0.05;    % This number times the stdev of baseline gamma power proxy 
%                         % plus the average baseline gamma power proxy 
%                         % will decide low vs. high gamma during LED
%                         % pulse                       

precedingInds=floor(precedingWindow(1)*LFP_Fs)+1:floor(precedingWindow(2)*LFP_Fs);
duringInds=floor(duringWindow(1)*LFP_Fs)+1:floor(duringWindow(2)*LFP_Fs);
baselineInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);

% Get gamma power proxy
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);

totalTrialLength=(1/LFP_Fs)*size(bandPassedLFPbySweep,2);
% Gamma power proxy stats.
disp('mean for entire trial');
disp(mean(mean(bandPassedLFPbySweep(:,:),2),1));
disp('stdev for entire trial');
disp(std(mean(bandPassedLFPbySweep(:,:),2)));
disp('LED ON average baseline');
ledONbaselineAv=mean(mean(bandPassedLFPbySweep(ledForSweeps>0,baselineInds),1),2);
disp(ledONbaselineAv);
disp('LED ON baseline stdev');
ledONbaselineStdev=std(mean(bandPassedLFPbySweep(ledForSweeps>0,baselineInds),2),1);
disp(ledONbaselineStdev);
disp('LED OFF average baseline');
ledOFFbaselineAv=mean(mean(bandPassedLFPbySweep(ledForSweeps==0,baselineInds),1),2);
disp(ledOFFbaselineAv);
disp('LED OFF baseline stdev');
ledOFFbaselineStdev=std(mean(bandPassedLFPbySweep(ledForSweeps==0,baselineInds),2),1);
disp(ledOFFbaselineStdev);

% Classify trials
trialClassifications=zeros(size(bandPassedLFPbySweep,1),1);
for i=1:size(bandPassedLFPbySweep,1)
    if ledForSweeps(i)==0
        if mean(bandPassedLFPbySweep(i,precedingInds))<=precedingGammaThresh*ledOFFbaselineStdev+ledOFFbaselineAv
            if mean(bandPassedLFPbySweep(i,duringInds))>duringGammaThresh*ledOFFbaselineStdev+ledOFFbaselineAv
                trialClassifications(i)=1;
            else
                trialClassifications(i)=2;
            end
        else
            if mean(bandPassedLFPbySweep(i,duringInds))>duringGammaThresh*ledOFFbaselineStdev+ledOFFbaselineAv
                trialClassifications(i)=3;
            else
                trialClassifications(i)=4;
            end
        end
    else % ledForSweeps(i)==1
        if mean(bandPassedLFPbySweep(i,precedingInds))<=precedingGammaThresh*ledONbaselineStdev+ledONbaselineAv
            if mean(bandPassedLFPbySweep(i,duringInds))>duringGammaThresh*ledONbaselineStdev+ledONbaselineAv
                trialClassifications(i)=5;
            else
                trialClassifications(i)=6;
            end
        else
            if mean(bandPassedLFPbySweep(i,duringInds))>duringGammaThresh*ledONbaselineStdev+ledONbaselineAv
                trialClassifications(i)=7;
            else
                trialClassifications(i)=8;
            end
        end
    end
end

% Display conditional probabilities for seeing high or low late activity
% Classification    LED      Gamma preceding        Gamma during
% 1                 off      low                    high
% 2                 off      low                    low
% 3                 off      high                   high
% 4                 off      high                   low
% 5                 on       low                    high
% 6                 on       low                    low
% 7                 on       high                   high
% 8                 on       high                   low
disp('Conditional probabilities');
disp('LED off: P(high late|low early)');
disp(sum(trialClassifications==1)/sum(trialClassifications==1|trialClassifications==2));
disp('LED off: P(low late|low early)');
disp(sum(trialClassifications==2)/sum(trialClassifications==1|trialClassifications==2));
disp('LED off: P(high late|high early)');
disp(sum(trialClassifications==3)/sum(trialClassifications==3|trialClassifications==4));
disp('LED off: P(low late|high early)');
disp(sum(trialClassifications==4)/sum(trialClassifications==3|trialClassifications==4));

disp('LED on: P(high late|low early)');
disp(sum(trialClassifications==5)/sum(trialClassifications==5|trialClassifications==6));
disp('LED on: P(low late|low early)');
disp(sum(trialClassifications==6)/sum(trialClassifications==5|trialClassifications==6));
disp('LED on: P(high late|high early)');
disp(sum(trialClassifications==7)/sum(trialClassifications==7|trialClassifications==8));
disp('LED on: P(low late|high early)');
disp(sum(trialClassifications==8)/sum(trialClassifications==7|trialClassifications==8));
        