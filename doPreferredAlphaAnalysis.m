function [matchOut,mismatchOut,cellPower]=doPreferredAlphaAnalysis(spikes,step1,s)

exptName='Mawake50';
duration=5;
binsize=1; % in ms
ledCondToDo=1; % 0 means both, 1 means LED off, 2 means LED on
saveDir='W:\Lea Final\Leas work\Lea\Getting Alpha Spikes\Alpha Preferred NonPreferred Analysis\No Base Subtract\';
togSaveDir=[saveDir exptName '_alphaPreffedAnalysis.mat'];
excludeTop4=1;
downSamp=30;

if excludeTop4==1
    maxinds=zeros(1,4);
    temp=step1.preffed.vals;
    for i=1:4
        [m,mi]=max(temp);
        temp(mi)=-1000;
        maxinds(i)=mi;
    end
    useInds=~ismember(1:length(step1.preffed.stim),maxinds);
    step1.PSTHs.ledOff=step1.PSTHs.ledOff(useInds,:,:);
    step1.PSTHs.ledOn=step1.PSTHs.ledOn(useInds,:,:);
    step1.PSTHs.ledBoth=step1.PSTHs.ledBoth(useInds,:,:);
    step1.preffed.stim=step1.preffed.stim(useInds);
    step1.preffed.vals=step1.preffed.vals(useInds);
    step1.nonpreffed.stim=step1.nonpreffed.stim(useInds);
    step1.nonpreffed.vals=step1.nonpreffed.vals(useInds);
end
    
if isempty(step1)
    s=unique(spikes.stimcond);
    f=unique(spikes.fileInd);
    [PSTHs,preffed,nonpreffed,x]=testAlpha(spikes,f,s,binsize,duration);
    step1.PSTHs=PSTHs;
    step1.preffed=preffed;
    step1.nonpreffed=nonpreffed;
    step1.x=x;
    save(togSaveDir,'step1');
else
    PSTHs=step1.PSTHs;
    preffed=step1.preffed;
    nonpreffed=step1.nonpreffed;
    x=step1.x;
end

[matchOut,mismatchOut,cellPower,cellByCell]=plotPreffedVsNonpreffedAlpha(PSTHs,preffed,nonpreffed,s,ledCondToDo,binsize);

figure(); 
plot(downSampAv(matchOut,downSamp),'Color','b');
hold on;
plot(downSampAv(mismatchOut,downSamp),'Color','g');

disp('done');