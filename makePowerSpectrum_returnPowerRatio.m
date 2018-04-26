function powerRatio=makePowerSpectrum_returnPowerRatio(data,Fs,lowerBand,upperBand)
% Trials are different rows

data=data';
% Now each column is a trial for fft

% Parameters
%Fs=6400;
T=1/Fs;
L=size(data,1);

nfft=2^nextpow2(L);

avSpec=[];
for i=1:size(data,2)
%     disp(i);
    Y=fft(data(:,i),nfft)/L;
    f=Fs/2*linspace(0,1,nfft/2+1);
    if i==1
        avSpec=2*abs(Y(1:nfft/2+1));
        allSpec(1,:)=2*abs(Y(1:nfft/2+1));
    else
        avSpec=avSpec+2*abs(Y(1:nfft/2+1));
        allSpec(i,:)=2*abs(Y(1:nfft/2+1));
    end
end
avSpec=avSpec/size(data,2);

% Plot average spectrum
figure;
plot(f,avSpec);
title('Averaged power spectrum');

powerRatio=sum(allSpec(:,f>=lowerBand(1) & f<lowerBand(2)),2)./sum(allSpec(:,f>=upperBand(1) & f<upperBand(2)),2);