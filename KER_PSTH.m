function KER_PSTH(spikes,trialLength,nParams1,param1Name,nParams2,param2Name,ONstart,ONstop,OFFstart,OFFstop)

[psth,x]=hist(spikes.spiketimes,0.05:0.1:trialLength);
plot(x,psth/(0.1*(nParams1*nParams2)),'Color','k');
title('Total PSTH');
xlabel('Time from Stimulus Onset (s)');
ylabel('F.R.(Hz)');

lineColors=zeros(100,3);
l=1;
for m=0:0.16:0.8
    for j=0:0.16:0.8
        for k=0:0.16:0.8
            lineColors(l,:)=[m j k];
            l=l+1;
        end
    end
end
lineColors=lineColors(randperm(size(lineColors,1)),:);
lineColorN=1;

h1=figure;
maxVal=1;
subPs=cell(nParams1,1);
% for i=1:nParams1
%     subPs{i}=subplot(nParams1,1,i);
%     for j=1:nParams2
%         hold on;
%         someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*nParams1+j);
%         [psth,x]=hist(someSpikes.spiketimes,0.25:0.5:trialLength);
%         plot(x,psth/0.5,'Color',lineColors(lineColorN,:));
%         lineColorN=lineColorN+1;
%         m=max(psth/0.5);
%         if m>maxVal
%             maxVal=m;
%         end
%     end 
% end
for i=1:nParams1
    subPs{i}=subplot(nParams1,1,i);
    for j=1:nParams2
        hold on;
        someSpikes=filtspikes(spikes,0,'stimcond',(i-1)*nParams1+j);
        [psth,x]=hist(someSpikes.spiketimes,0.1:0.2:trialLength);
        plot(x,psth/0.2,'Color',lineColors(lineColorN,:));
        lineColorN=lineColorN+1;
        %plot([2.5 2.5],[0 30],'Color','k');
        m=max(psth/0.2);
        if m>maxVal
            maxVal=m;
        end
    end 
end

for i=1:nParams1
    set(h1,'CurrentAxes',subPs{i});
    axis([0 trialLength 0 maxVal]);
    if i==1
        title(['PSTH - Different Subplots: ' param1Name ', Different Lines: ' param2Name]);
        xlabel('Time from Stimulus Onset (s)');
        ylabel('F.R.(Hz)');
    end
end

h2=figure;
maxVal=1;
clear subPs
subPs=cell(nParams1,1);
% for i=1:nParams1
%     subPs{i}=subplot(nParams1,1,i);
%     hold on;
%     someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+nParams2-1);
%     [psth,x]=hist(someSpikes.spiketimes,0.025:0.05:trialLength);
%     plot(x,psth/(0.05*nParams2),'Color',lineColors(lineColorN,:));
%     lineColorN=lineColorN+1;
%     m=max(psth/(0.05*nParams2));
%     if m>maxVal
%         maxVal=m;
%     end
% end
for i=1:nParams1
    subPs{i}=subplot(nParams1,1,i);
    hold on;
    someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+nParams2-1);
    [psth,x]=hist(someSpikes.spiketimes,0.01:0.02:trialLength);
    plot(x,psth/(0.02*nParams2),'Color',lineColors(lineColorN,:));
    lineColorN=lineColorN+1;
    m=max(psth/(0.02*nParams2));
    if m>maxVal
        maxVal=m;
    end
end
for i=1:nParams1
    set(h2,'CurrentAxes',subPs{i});
    axis([0 trialLength 0 maxVal]);
    if i==1
        title(['PSTH - Different Subplots For Diff. Values of the Variable: ' param1Name]);
        xlabel('Time from Stimulus Onset (s)');
        ylabel('F.R. (Hz)');
    end
end

figure;
lineColorN=1;
for i=1:nParams1
    someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+nParams2-1);
    [psth,x]=hist(someSpikes.spiketimes,0.05:0.1:trialLength);
    plot(x,psth/(0.1*nParams2),'Color',lineColors(lineColorN,:));
    lineColorN=lineColorN+1;
    m=max(psth/(0.1*nParams2));
    if m>maxVal
        maxVal=m;
    end
    hold on;
end
%axis([0 trialLength 0 maxVal]);
title(['PSTH - Different Subplots For Diff. Values of the Variable: ' param1Name]);
xlabel('Time from Stimulus Onset (s)');
ylabel('F.R. (Hz)');

figure;
ONresponse=zeros(nParams1,1);
OFFresponse=zeros(nParams1,1);
for i=1:nParams1
    subplot(nParams1,1,i);
    hold on;
    someSpikes=filtspikes(spikes,0,'stimcond',i:1:i+nParams2-1);
    [psth2,x2]=hist(someSpikes.spiketimes,0.05:0.1:trialLength);
    ONresponse(i)=max(psth2((x2>ONstart) & (x2<ONstop)));
    OFFresponse(i)=max(psth2((x2>OFFstart) & (x2<OFFstop)));
end

if strcmp(param1Name,'Orientation')
    figure;
    orientStep=360/nParams1;
    orients=0:orientStep:360;
    ONresponse=[ONresponse; ONresponse(1)];
    OFFresponse=[OFFresponse; OFFresponse(1)];
    ONnormMags=ONresponse/max(ONresponse);
    OFFnormMags=OFFresponse/max(OFFresponse);
    subplot(1,2,1);
    if max(ONresponse)>max(OFFresponse)
        polar(deg2rad(orients)',ONresponse,'-g');
        hold on;
        title('Orientation Tuning');
        polar(deg2rad(orients)',OFFresponse,'-r');
        subplot(1,2,2);
        polar(deg2rad(orients)',ONnormMags,'-g');
        hold on;
        polar(deg2rad(orients)',OFFnormMags,'-r');
        title('ON: Green, OFF: Red');
    else
        polar(deg2rad(orients)',OFFresponse,'-r');
        hold on;
        title('Orientation Tuning');
        polar(deg2rad(orients)',ONresponse,'-g');
        subplot(1,2,2);
        polar(deg2rad(orients)',OFFnormMags,'-r');
        hold on;
        polar(deg2rad(orients)',ONnormMags,'-g');
        title('ON: Green, OFF: Red');
    end
else
    figure;
    bar(1:nParams1,ONresponse');
    title('ON Response for Different Values of Stimulus Condition (Variable 1)');
    figure;
    bar(1:nParams1,OFFresponse');
    title('OFF Response for Different Values of Stimulus Condition (Variable 1)');
end

end


