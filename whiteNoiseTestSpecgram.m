function whiteNoiseTestSpecgram(Fs,numSamples)

x=0:1/Fs:(1/Fs)*(numSamples-1);
freqs=[0.1:0.1:0.9 1:100];
y=zeros(1,length(x));
for i=1:length(freqs)
    y=y+sin(2*pi*freqs(i).*x);
    y=y+sin(2*pi*freqs(i).*x+(2*pi*(1/6)));
    y=y+sin(2*pi*freqs(i).*x+(2*pi*(2/6)));
    y=y+sin(2*pi*freqs(i).*x+(2*pi*(3/6)));
    y=y+sin(2*pi*freqs(i).*x+(2*pi*(4/6)));
    y=y+sin(2*pi*freqs(i).*x+(2*pi*(5/6)));
end

[p,LFPspecgram]=makeWaveletSpecgram(y,1,Fs);
