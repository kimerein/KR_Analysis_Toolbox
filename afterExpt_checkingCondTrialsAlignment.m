% Color LED data by LEDConds to make sure trials are aligned
rigdef=RigDefs();

figure();
disp('Making LED figure');
[Fs,ledAv,ledData,gotDaqs,allCompleteSweeps,fromFileInds]=getLEDbySweep(rigdef.Dir.Data,expt.files.names,rigdef.Daq.SampleRate,10,rigdef.Daq.LEDChannel);
j=1;
for i=1:length(expt.sweeps.led)
%     if a(j)~=expt.sweeps.fileInd(i)
%         continue;
%     end
    if expt.sweeps.led(i)==0
        plot(ledData(j,:),'Color','black');
    elseif expt.sweeps.led(i)>0 && expt.sweeps.led(i)<2
        plot(ledData(j,:),'Color','blue');
    elseif expt.sweeps.led(i)>=2 && expt.sweeps.led(i)<3
        plot(ledData(j,:),'Color','cyan');
    elseif expt.sweeps.led(i)>=3 && expt.sweeps.led(i)<4
        plot(ledData(j,:),'Color','gray');
    elseif expt.sweeps.led(i)>=4 && expt.sweeps.led(i)<=5
        plot(ledData(j,:),'Color','red');
    end
    hold on;
    j=j+1;
end