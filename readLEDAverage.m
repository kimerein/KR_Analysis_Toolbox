function ledAv=readLEDAverage(daqFileName,Fs,ledChannel)

global dataDir

data=daqread(strcat(dataDir,'\',daqFileName));
led