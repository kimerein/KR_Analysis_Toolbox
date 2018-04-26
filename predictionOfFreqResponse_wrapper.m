function [freq,amps]=predictionOfFreqResponse_wrapper(timecourse,timecourse_anesth,thal_freq,thalamusResponse,transFunc)
% make sure input x is in ms, not seconds
% beware photoartifact
% downSampFactor=1;
% timecourse.x=downSampAv(timecourse.x,downSampFactor);
% timecourse.y=downSampAv(timecourse.y,downSampFactor);
% timecourse_anesth.x=downSampAv(timecourse_anesth.x,downSampFactor);
% timecourse_anesth.y=downSampAv(timecourse_anesth.y,downSampFactor);

% afterTauBaseline=[80 120]; % in ms
afterTauBaseline=[80 100]; % in ms
% afterTauBaseline=[80 100]; % in ms
% afterTauBaseline=[100 110]; % in ms
w=3;
% afterTauBaseline=[120 130]; % in ms
takeTauUntil=100;
% takeTauUntil=100;
% takeTauUntil=110;
dontBaseShiftFromAnesth=0;
putAtZero=0;
scaleToOne=0;
useOnlyAnesth=0;
putAtZeroAfterAll=0;

figure(); plot(timecourse_anesth.x,timecourse_anesth.y);
figure(); plot(timecourse.x,timecourse.y);

anesth_baseline=mean(timecourse_anesth.y(timecourse_anesth.x>=-500 & timecourse_anesth.x<=-300),2);
anesth_afterTau=mean(timecourse_anesth.y(timecourse_anesth.x>=afterTauBaseline(1) & timecourse_anesth.x<=afterTauBaseline(2)),2);
anesth_ev=mean(timecourse_anesth.y(timecourse_anesth.x>=-20 & timecourse_anesth.x<=-4),2);

all_baseline=mean(timecourse.y(timecourse.x>=-500 & timecourse.x<=-300),2);
all_afterTau=mean(timecourse.y(timecourse.x>=afterTauBaseline(1) & timecourse.x<=afterTauBaseline(2)),2);
all_ev=mean(timecourse.y(timecourse.x>=-20 & timecourse.x<=-4),2);

dip=(anesth_baseline-anesth_afterTau)/(anesth_ev-anesth_baseline);
allEvSize=all_ev-all_baseline;
allEvDip=allEvSize*dip;
subtractOff=allEvDip+all_afterTau;


useTimecourse=timecourse.y(timecourse.x>=w & timecourse.x<=takeTauUntil);
if putAtZero==1
    useTimecourse=useTimecourse-all_afterTau;
elseif dontBaseShiftFromAnesth==0
    useTimecourse=useTimecourse-subtractOff;
end
if useOnlyAnesth==1
    useTimecourse=timecourse_anesth.y(timecourse_anesth.x>=w & timecourse_anesth.x<=takeTauUntil);
end
if putAtZeroAfterAll==1
    useTimecourse=useTimecourse-mean(useTimecourse(end-5:end));
end

figure(); plot((timecourse.x(timecourse.x>=w & timecourse.x<=takeTauUntil)-w)./1000,useTimecourse);

[giveOutput,freq]=predictTFcutoff_withThalResponse((timecourse.x(timecourse.x>=w & timecourse.x<=takeTauUntil)-w)./1000,useTimecourse,[thal_freq(1:6) interp(thal_freq(7:end),1)],[mean(thalamusResponse(:,1:6),1) interp(mean(thalamusResponse(:,7:end),1),1)],transFunc);
% [giveOutput,freq]=predictTFcutoff_withThalResponse((timecourse.x(timecourse.x>=w & timecourse.x<=takeTauUntil)-w)./1000,useTimecourse,[thal_freq(1:6) interp(thal_freq(7:end),10)],[mean(thalamusResponse(:,1:6),1) interp(mean(thalamusResponse(:,7:end),1),10)],transFunc);
% [giveOutput,freq]=predictTFcutoff_withThalResponse((timecourse.x(timecourse.x>=w & timecourse.x<=takeTauUntil)-w)./1000,useTimecourse,interp(thal_freq,100),interp(mean(thalamusResponse,1),100),transFunc);
figure(); 
amps=sqrt(giveOutput);
semilogx(freq,amps); 


