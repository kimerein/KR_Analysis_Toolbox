function [matchOut,mismatchOut,cellPower,cellByCell]=plotPreffedVsNonpreffedAlpha(PSTHs,preffed,nonpreffed,useStimcondInd,ledCondToDo,binsize)
% binsize is in ms
timeWindow=[4 5]; % in sec relative to stimulus onset


switch ledCondToDo
    case 0
    case 1
        PSTHs.ledBoth=PSTHs.ledOff;
    case 2
        PSTHs.ledBoth=PSTHs.ledOn;
end

match=[];
mismatch=[];
for i=1:length(useStimcondInd)
    if i==1
        match=PSTHs.ledBoth(useStimcondInd(preffed.stim)==useStimcondInd(i),i,:);
        mismatch=PSTHs.ledBoth(useStimcondInd(nonpreffed.stim)==useStimcondInd(i),i,:);
    else
        match=[match; PSTHs.ledBoth(useStimcondInd(preffed.stim)==useStimcondInd(i),i,:)];
        mismatch=[mismatch; PSTHs.ledBoth(useStimcondInd(nonpreffed.stim)==useStimcondInd(i),i,:)];
    end
end

% Looking at individual cells
for i=1:size(PSTHs.ledBoth,1)
    % Go cell by cell
    cellmatch(i,:,:)=PSTHs.ledBoth(i,useStimcondInd==useStimcondInd(preffed.stim(i)),:);
    cellmismatch(i,:,:)=PSTHs.ledBoth(i,useStimcondInd==useStimcondInd(nonpreffed.stim(i)),:);
end
cellpowerForMatches=nan(size(cellmatch,1),1);
cellpowerForNonMatches=nan(size(cellmismatch,1),1);
cellpowerForMatches8to12=nan(size(cellmatch,1),1);
cellpowerForNonMatches8to12=nan(size(cellmismatch,1),1);
cellxtimes=linspace(0,5,size(cellmismatch,3));
for i=1:size(cellmatch,1)
    L=length(cellmatch(i,:,:));
    Fs=1/(binsize/1000);
    NFFT=2^nextpow2(L);
    useTimeInds=cellxtimes>=timeWindow(1) & cellxtimes<=timeWindow(2);
    Y=fft(cellmatch(i,useTimeInds),NFFT)/L;
    fff=Fs/2*linspace(0,1,NFFT/2+1);
    amp=2*abs(Y(1:NFFT/2+1));
%     normPower=mean(amp(fff>=0.5 & fff<=4));
    normPower=sum(amp);
%     normPower=1;
%     normPower=mean(amp(fff>=30 & fff<=100));
    if normPower==0
        continue
    end
    alphaPow=mean(amp(fff>=4 & fff<=8));
    cellpowerForMatches(i)=alphaPow/normPower;
    highAlphaPow=mean(amp(fff>=8 & fff<=12));
    cellpowerForMatches8to12(i)=highAlphaPow/normPower;
end
for i=1:size(cellmismatch,1)
    L=length(cellmismatch(i,:,:));
    Fs=1/(binsize/1000);
    NFFT=2^nextpow2(L);
    useTimeInds=cellxtimes>=timeWindow(1) & cellxtimes<=timeWindow(2);
    Y=fft(cellmismatch(i,useTimeInds),NFFT)/L;
    fff=Fs/2*linspace(0,1,NFFT/2+1);
    amp=2*abs(Y(1:NFFT/2+1));
%     normPower=mean(amp(fff>=0.5 & fff<=4));
    normPower=sum(amp);
%     normPower=1;
%     normPower=mean(amp(fff>=30 & fff<=100));
    if normPower==0
        continue
    end
    alphaPow=mean(amp(fff>=4 & fff<=8));
    cellpowerForNonMatches(i)=alphaPow/normPower;
    highAlphaPow=mean(amp(fff>=8 & fff<=12));
    cellpowerForNonMatches8to12(i)=highAlphaPow/normPower;
end
cellByCell.matches4to8=cellpowerForMatches;
cellByCell.nonmatches4to8=cellpowerForNonMatches;
cellByCell.matches8to12=cellpowerForMatches8to12;
cellByCell.nonmatches8to12=cellpowerForNonMatches8to12;    

% for i=1:length(useStimcondInd)
%     if i==1
%         match=PSTHs.ledOn(preffed.stim==useStimcondInd(i),useStimcondInd(i),:);
%         mismatch=PSTHs.ledOn(nonpreffed.stim==useStimcondInd(i),useStimcondInd(i),:);
%     else
%         match=[match; PSTHs.ledOn(preffed.stim==useStimcondInd(i),useStimcondInd(i),:)];
%         mismatch=[mismatch; PSTHs.ledOn(nonpreffed.stim==useStimcondInd(i),useStimcondInd(i),:)];
%     end
% end

% match=nanmean(match);
% mismatch=nanmean(mismatch);
powerForMatches=nan(size(match,1),1);
powerForNonMatches=nan(size(mismatch,1),1);
powerForMatches8to12=nan(size(match,1),1);
powerForNonMatches8to12=nan(size(mismatch,1),1);
xtimes=linspace(0,5,size(mismatch,3));
for i=1:size(match,1)
    L=length(match(i,:,:));
    Fs=1/(binsize/1000);
    NFFT=2^nextpow2(L);
    useTimeInds=find(xtimes>=timeWindow(1) & xtimes<=timeWindow(2));
    Y=fft(reshape(match(i,:,useTimeInds),1,length(useTimeInds)),NFFT)/L;
    fff=Fs/2*linspace(0,1,NFFT/2+1);
    amp=2*abs(Y(1:NFFT/2+1));
%     normPower=mean(amp(fff>=0.5 & fff<=4));
    normPower=sum(amp);
%     normPower=1;
%     normPower=mean(amp(fff>=30 & fff<=100));
    if normPower==0
        continue
    end
    alphaPow=mean(amp(fff>=4 & fff<=8));
    powerForMatches(i)=alphaPow/normPower;
    highAlphaPow=mean(amp(fff>=8 & fff<=12));
    powerForMatches8to12(i)=highAlphaPow/normPower;
end
for i=1:size(mismatch,1)
    L=length(mismatch(i,:,:));
    Fs=1/(binsize/1000);
    NFFT=2^nextpow2(L);
    useTimeInds=find(xtimes>=timeWindow(1) & xtimes<=timeWindow(2));
    Y=fft(reshape(mismatch(i,:,useTimeInds),1,length(useTimeInds)),NFFT)/L;
    fff=Fs/2*linspace(0,1,NFFT/2+1);
    amp=2*abs(Y(1:NFFT/2+1));
%     normPower=mean(amp(fff>=0.5 & fff<=4));
    normPower=sum(amp);
%     normPower=1;
%     normPower=mean(amp(fff>=30 & fff<=100));
    if normPower==0
        continue
    end
    alphaPow=mean(amp(fff>=4 & fff<=8));
    powerForNonMatches(i)=alphaPow/normPower;
    highAlphaPow=mean(amp(fff>=8 & fff<=12));
    powerForNonMatches8to12(i)=highAlphaPow/normPower;
end
cellPower.matches4to8=powerForMatches;
cellPower.nonmatches4to8=powerForNonMatches;
cellPower.matches8to12=powerForMatches8to12;
cellPower.nonmatches8to12=powerForNonMatches8to12;    


matchOut=reshape(nanmean(match),1,size(match,3));
mismatchOut=reshape(nanmean(mismatch),1,size(mismatch,3));
% matchOut=reshape(match,1,size(match,3));
% mismatchOut=reshape(mismatch,1,size(mismatch,3));
figure(); 
plot(downSampAv(matchOut,10),'Color','b');
hold on;
plot(downSampAv(mismatchOut,10),'Color','g');