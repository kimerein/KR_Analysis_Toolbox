function [histPlots,lambdas]=fitPoissonToFreqResponse(freqs,freqResponse)

scaleFactor=2*287;
nBins=20;
% binStep=0.025;
binStep=0.025;
outliercutoff=10*binStep; % V1
% outliercutoff=10*binStep; % LGN
% outliercutoff=40*binStep; % prediction
poissFactor=1/binStep;
% binCenters=-0.5:binStep:1; % LGN
binCenters=-0.5:binStep:2; % V1

freqResponse=freqResponse./scaleFactor;

lambdas=nan(1,length(freqs));
for i=1:length(freqs)
    [n,xout]=hist(freqResponse(:,i),binCenters);
    histPlots{i}.xout=xout.*poissFactor;
    histPlots{i}.n=n;
    temp=freqResponse(:,i);
    temp(temp<0)=0;
    if ~isempty(outliercutoff)
        temp(temp>outliercutoff)=nan;
        temp=temp(~isnan(temp));
    end
    lambdahat=poissfit(temp.*poissFactor);
    histPlots{i}.temp=temp.*poissFactor;
    histPlots{i}.freqResponse=freqResponse(:,i).*poissFactor;
    lambdas(i)=lambdahat;
end

figure();
for i=1:length(freqs)
    subplot(15,1,i);
    plot(histPlots{i}.xout,histPlots{i}.n);
    hold on;
    y=poisspdf(floor(histPlots{i}.xout),lambdas(i));
    y=(y./max(y)).*max(histPlots{i}.n);
    plot(histPlots{i}.xout,y,'Color','r');
end

figure(); 
semilogx(freqs,lambdas);
    