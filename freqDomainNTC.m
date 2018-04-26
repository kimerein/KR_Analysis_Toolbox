function fromNTC=freqDomainNTC(timeDomain_x,timeDomain_y,expNTC)
% timeDomain_x must be in ms

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

if ~isempty(expNTC)
    timeDomain_x=linspace(0,2,2000);
    timeDomain_y=exp(-timeDomain_x./expNTC);
end

% Zero pad
timeDomain_x=linspace(0,4,4000);
timeDomain_y=[timeDomain_y zeros(1,4000-length(timeDomain_y))];

% Transform to frequency domain
L=length(timeDomain_x);
Fs=1/0.001; % in ms
NFFT=2^nextpow2(L);
Y=fft(timeDomain_y,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum
figure(); 
semilogx(f,2*abs(Y(1:NFFT/2+1)));

% Get frequency domain curve
useF=[];
useFInd=[];
for i=1:length(freqs)
    useFInd(i)=find(f>=freqs(i),1,'first');
    useF(i)=f(useFInd(i));
end

% Plot frequency domain curve
figure(); 
semilogx(f(useFInd),2*Y(useFInd));

% Store results
fromNTC.freqs=f(useFInd);
fromNTC.F_NTC=2*abs(Y(useFInd));
fromNTC.full_freqs=f;
fromNTC.full_F_NTC=2*abs(Y(1:NFFT/2+1));
fromNTC.timeDomain_x=timeDomain_x;
fromNTC.timeDomain_y=timeDomain_y;
