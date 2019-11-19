function spikes=fixSpikesAfterO2Sort(spikes,expt,detectfilenames)

sweepsPerDaqFile=4;

spikes.sweeps.fromExpt_fileInd=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));
spikes.sweeps.trigger=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));
spikes.sweeps.fileInd=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));
spikes.sweeps.trials=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));
spikes.sweeps.stimcond=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));
spikes.sweeps.led=single(ones(1,length(detectfilenames)*sweepsPerDaqFile));

for i=1:length(detectfilenames)
    for j=1:length(expt.files.names)
        if strcmp(expt.files.names{j},detectfilenames{i})
            spikes.sweeps.fromExpt_fileInd((i-1)*sweepsPerDaqFile+1:(i-1)*sweepsPerDaqFile+sweepsPerDaqFile)=j;
            spikes.sweeps.trigger((i-1)*sweepsPerDaqFile+1:(i-1)*sweepsPerDaqFile+sweepsPerDaqFile)=1:sweepsPerDaqFile;
        end
    end
end

spikes.sweeps.fileInd=spikes.sweeps.fromExpt_fileInd-nanmin(spikes.sweeps.fromExpt_fileInd)+1;
spikes.sweeps.trials=unique(spikes.trials);

spikes.led=single(ones(size(spikes.trials)));
spikes.stimcond=single(ones(size(spikes.trials)));
spikes.fileInd=single(ones(size(spikes.trials)));

for i=1:length(spikes.fileInd)
    spikes.fileInd(i)=spikes.sweeps.fileInd(spikes.sweeps.trials==spikes.trials(i));
    spikes.fromExpt_fileInd(i)=spikes.sweeps.fromExpt_fileInd(spikes.sweeps.trials==spikes.trials(i));
    spikes.trigger(i)=spikes.sweeps.trigger(spikes.sweeps.trials==spikes.trials(i));
end

spikes=addFieldtoSpikes(spikes,expt,'led');
spikes=addFieldtoSpikes(spikes,expt,'stimcond');

end

function newSpikes=addFieldtoSpikes(spikes,expt,fieldname)

newSpikes=spikes;
tempFieldVals=newSpikes.(fieldname);
exptFieldVals=expt.sweeps.(fieldname);
for i=1:length(spikes.(fieldname))
    currFileInd=spikes.fromExpt_fileInd(i);
    currTrigger=spikes.trigger(i);
    % Find condition for this daq file for this trigger
    firstIndForThisDaqFile=find(expt.sweeps.fileInd==currFileInd,1,'first');
    corrIndIntoLED=firstIndForThisDaqFile+(currTrigger-1);
    tempFieldVals(i)=exptFieldVals(corrIndIntoLED);
end
newSpikes.(fieldname)=tempFieldVals;

% Now fix spikes.sweeps
tempFieldVals=newSpikes.sweeps.(fieldname);
for i=1:length(newSpikes.sweeps.led)
    currFileInd=newSpikes.sweeps.fromExpt_fileInd(i);
    currTrigger=newSpikes.sweeps.trigger(i);
    firstInd=find(expt.sweeps.fileInd==currFileInd,1,'first');
    corrIndIntoLED=firstInd+(currTrigger-1);
    tempFieldVals(i)=exptFieldVals(corrIndIntoLED);
end
newSpikes.sweeps.(fieldname)=tempFieldVals;
   
end

