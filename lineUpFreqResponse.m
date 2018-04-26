function lineUpFreqResponse(freq1,resp1,freq2,resp2)

if length(freq1)>length(freq2)
    tempfreq=freq1;
    tempresp=resp1;
    freq1=freq2;
    resp1=resp2;
    freq2=tempfreq;
    resp2=tempresp;
end
mean1=mean

tryGains=1:1