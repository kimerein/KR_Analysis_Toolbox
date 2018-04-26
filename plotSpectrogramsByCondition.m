function fig_handles=plotSpectrogramsByCondition(LFPdata,LFP_Fs,stimsForSweeps,ledForSweeps,params,matchNumTrials,useNTrials,makeSpecSets)
% KR - make spectrogram with wavelets (code from Ed's lab rotation)
% useNTrials gives the max number of trials to use to make each averaged
% spectrogram
% Increasing useNTrials substantially increases running time
% makeSpecSets(1)=1 if want to make spectrogram set 1 (LED ON vs. OFF,
% collapsed over all stim. conditions); else makeSpecSets(1)=0
% makeSpecSet(2)=1 if want to make spectrogram set 2 (LED ON vs. OFF for
% specific stimulus conditions)

fig_handles=[];
newF=[];

% Normalize ledForSweeps to make it a logical array
% Assumes a step LED (that is, only 2 LED values)
minLED=min(ledForSweeps);
maxLED=max(ledForSweeps);
if minLED==maxLED
    ledForSweeps(ledForSweeps==0)=0;
    ledForSweeps(ledForSweeps~=0)=1;
else
    ledForSweeps(ledForSweeps==minLED)=0;
    ledForSweeps(ledForSweeps==maxLED)=1;
end


% SPECTROGRAM SET 1
disp('Making spectrograms');
set1=find(logical(ledForSweeps));
set2=find(~logical(ledForSweeps));
if matchNumTrials==1
    [set1,set2]=matchNumberOfTrials(set1,set2);
end

if ~isnumeric(useNTrials) % 'all'
else
    if useNTrials>length(set1)
        useN=length(set1);
    else
        useN=useNTrials;
    end
    these=randperm(length(set1));
    set1=set1(these(1:useN));
    if useNTrials>length(set2)
        useN=length(set2);
    else
        useN=useNTrials;
    end
    these=randperm(length(set2));
    set2=set2(these(1:useN));
end

if makeSpecSets(1)==1 && ~isempty(set1) && ~isempty(set2)
    if any(set1~=0)
        set1
        newF=makeWaveletSpecgram(LFPdata,set1,LFP_Fs);
    end
    if ~isempty(newF)
        title('Spectrogram - Collapsed Over All Stim Conds. - LED ON');
        fig_handles=[fig_handles; newF];
    end
    if any(set2~=0)
        set2
        newF=makeWaveletSpecgram(LFPdata,set2,LFP_Fs);
    end
    if ~isempty(newF)
        title('Spectrogram - Collapsed Over All Stim Conds. - LED OFF');
        fig_handles=[fig_handles; newF];
    end
end


% SPECTROGRAM SET 2
if makeSpecSets(2)==1
    for i=1:length(params.Var2_values)
        disp('Making spectrograms');
        theseConds=i:length(params.Var2_values):length(params.Var1_values)*length(params.Var2_values);
        set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
        set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
        if matchNumTrials==1
            [set1,set2]=matchNumberOfTrials(set1,set2);
        end
        if ~isnumeric(useNTrials) % 'all'
        else
            if useNTrials>length(set1)
                useN=length(set1);
            else
                useN=useNTrials;
            end
            these=randperm(length(set1));
            set1=set1(these(1:useN));
            if useNTrials>length(set2)
                useN=length(set2);
            else
                useN=useNTrials;
            end
            these=randperm(length(set2));
            set2=set2(these(1:useN));
        end  
        if isempty(set1) || isempty(set2)
            continue;
        end
        if any(set1~=0)
            newF=makeWaveletSpecgram(LFPdata,set1,LFP_Fs);
        end
        if ~isempty(newF)
            title(['Spectrogram - Collapsed Over ' params.Var1_name ' - LED ON']);
            fig_handles=[fig_handles; newF];
        end
        if any(set2~=0)
            newF=makeWaveletSpecgram(LFPdata,set2,LFP_Fs);
        end
        if ~isempty(newF)
            title(['Spectrogram - Collapsed Over ' params.Var1_name ' - LED OFF']);
            fig_handles=[fig_handles; newF];
        end
    end
end
if makeSpecSets(2)==1
    for i=1:length(params.Var1_values)
        disp('Making spectrograms');
        theseConds=i:1:i+length(params.Var2_values)-1;
        set1=ismember(stimsForSweeps,theseConds)&logical(ledForSweeps);
        set2=ismember(stimsForSweeps,theseConds)&~logical(ledForSweeps);
        if matchNumTrials==1
            [set1,set2]=matchNumberOfTrials(set1,set2);
        end
        if ~isnumeric(useNTrials) % 'all'
        else
            if useNTrials>length(set1)
                useN=length(set1);
            else
                useN=useNTrials;
            end
            these=randperm(length(set1));
            set1=set1(these(1:useN));
            if useNTrials>length(set2)
                useN=length(set2);
            else
                useN=useNTrials;
            end
            these=randperm(length(set2));
            set2=set2(these(1:useN));
        end  
        if isempty(set1) || isempty(set2)
            continue;
        end
        if any(set1~=0)
            newF=makeWaveletSpecgram(LFPdata,set1,LFP_Fs);
        end
        if ~isempty(newF)
            title(['Spectrogram - Collapsed Over ' params.Var2_name ' - LED ON']);
            fig_handles=[fig_handles; newF];
        end
        if any(set2~=0)
            newF=makeWaveletSpecgram(LFPdata,set2,LFP_Fs);
        end
        if ~isempty(newF)
            title(['Spectrogram - Collapsed Over ' params.Var2_name ' - LED OFF']);
            fig_handles=[fig_handles; newF];
        end
    end
end       

function [set1,set2]=matchNumberOfTrials(set1,set2)
    if length(set1)<length(set2)
        set2=set2(randi(length(set2),1,length(set1)));
    elseif length(set2)<length(set1)
        set1=set1(randi(length(set1),1,length(set2)));
    end
