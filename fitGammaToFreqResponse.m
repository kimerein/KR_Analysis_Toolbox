function fitGammaToFreqResponse(freqs,amps)

modelFun=@(k,theta) (x.^(k-1).*exp(-x./theta))./(theta.^k.*gamma(k));
startingVals=