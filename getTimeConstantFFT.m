function fromNTC=getTimeConstantFFT(x,y,zeroSubtract)

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

if zeroSubtract==1
    y=y-mean(y(end-10:end));
end

% figure(); 
% plot(x,y);

% x=linspace(0,2,10000);
% y=[y zeros(1,10000-length(y))];
x=linspace(0,10,10000);
y=[y zeros(1,10000-length(y))];

Fs=1/(x(2)-x(1));
L = length(x);        % Length of signal

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
% figure();
% semilogx(f,2*abs(Y(1:NFFT/2+1))) ;
% title('Single-Sided Amplitude Spectrum of y(t)');
% xlabel('Frequency (Hz)');
% ylabel('|Y(f)|');

useF=[]; useFInd=[];
for i=1:15
useFInd(i)=find(f>=freqs(i),1,'first');
useF(i)=f(useFInd(i));
end

fromNTC.full_F_NTC=2*abs(Y(1:NFFT/2+1));
fromNTC.freqs=f(useFInd);
fromNTC.F_NTC=2*abs(Y(useFInd));
fromNTC.full_freqs=f;
