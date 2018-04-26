function [ssAmps,normSsAmps]=getThalFieldsAcrossFreqs(x,L,ledConds)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% ledVals=freqs+0.05;
ledVals=freqs;

ssAmps=zeros(1,length(freqs));
for i=1:length(freqs)
    currFreq=freqs(i);
    [amps,~,normAmps]=alignTracesAndPlotAverage(x,mean(L(ismember(ledConds,ledVals(i)),:),1),[],freqs(i));
%     if length(amps)>2
%         ssAmps(i)=mean(amps(3:end));
%     else
%         ssAmps(i)=amps(2);
%     end
%     if length(normAmps)>2
%         normSsAmps(i)=mean(normAmps(3:end));
%     else
%         normSsAmps(i)=normAmps(2);
%     end
    if length(amps)>5
        ssAmps(i)=mean(amps(6:end));
    else
        ssAmps(i)=mean(amps(2:end));
    end
    if length(normAmps)>5
        normSsAmps(i)=mean(normAmps(6:end));
    else
        normSsAmps(i)=mean(normAmps(2:end));
    end
end

figure();
plot(freqs,ssAmps);

figure();
plot(freqs,normSsAmps);