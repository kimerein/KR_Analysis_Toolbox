function makePSTHalignedToLED(psth,LEDbySweep,ledConds,ledVal,trialDuration,saveDir)

if iscell(LEDbySweep)
    LEDbySweep=LEDbySweep{1};
end

if size(LEDbySweep,1)~=size(psth.psths{1},1)
    error('Sizes of LEDbySweep and psth do not match');
end
if ~isequal(ledConds,psth.unitLED{1})
    error('LED conditions for LEDbySweep and psth do not match');
end

% Align psth to led
times=linspace(0,trialDuration,size(LEDbySweep,2));
% [psth,ledConds,ledStart]=alignPSTHtoLED(psth,LEDbySweep,ledConds,ledVal,times,[]);
[psth,ledConds,ledStart]=alignPSTHtoLED(psth,LEDbySweep,ledConds,ledVal,times,3);

save([saveDir '\psth_alignedToLED.mat'],'psth');
save([saveDir '\ledConds_alignedToLED.mat'],'ledConds');
save([saveDir '\ledStart_alignedToLED.mat'],'ledStart');