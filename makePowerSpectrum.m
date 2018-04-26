function [f,avSpec,allSpecs]=makePowerSpectrum(data,Fs)

showFigs=0;
% Normalize or whiten?
normalize=0;
normBand=[1 70]; % in Hz
whiten=0;

% Trials are different rows
data=data';
% Now each column is a trial for fft

% Parameters
%Fs=6400;
T=1/Fs;
L=size(data,1);

nfft=2^nextpow2(L);

avSpec=[];
allSpecs=[];
for i=1:size(data,2)
%     disp(i);
    Y=fft(data(:,i),nfft)/L;
    f=Fs/2*linspace(0,1,nfft/2+1);
    if i==1
        avSpec=2*abs(Y(1:nfft/2+1));
        allSpecs(i,:)=2*abs(Y(1:nfft/2+1));
    else
        avSpec=avSpec+2*abs(Y(1:nfft/2+1));
        allSpecs(i,:)=2*abs(Y(1:nfft/2+1));
    end
end
avSpec=avSpec/size(data,2);

% Plot average spectrum
if showFigs==1
    figure;
    plot(f,avSpec);
    xlim(normBand);
    title('Averaged power spectrum');
end

if normalize==1
    figure(); 
    plot(f(f>=normBand(1) & f<=normBand(2)),avSpec(f>=normBand(1) & f<=normBand(2))./sum(avSpec(f>=normBand(1) & f<=normBand(2))));  
end

if whiten==1
    if normalize~=1
        avSpec=avSpec.*((f').^1);
        if showFigs==1
            figure(); 
            plot(f,avSpec.*((f').^1));
        end
    else
        figure(); 
        temp=avSpec.*((f').^1);
%         temp=(avSpec(f>=normBand(1) & f<=normBand(2)))./(f(f>=normBand(1) & f<=normBand(2)))';
        plot(f(f>=normBand(1) & f<=normBand(2)),temp(f>=normBand(1) & f<=normBand(2))./sum(temp(f>=normBand(1) & f<=normBand(2))));
    end
end