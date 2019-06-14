function newspikes=prep_filtspikes(spikes,code,ledValue,stimValue)

switch code
    case 'led'
        temp=[];
        temp1=[];
        for i=1:length(ledValue)
            spikes=makeTempField(spikes,'led',ledValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        newspikes=filtspikes(spikes,0,'temp',1);
    case 'stimcond'
        temp=[];
        temp1=[];
        for i=1:length(stimValue)
            spikes=makeTempField(spikes,'stimcond',stimValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        newspikes=filtspikes(spikes,0,'temp',1);
end