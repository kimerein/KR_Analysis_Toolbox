function bandPassedLFP=bandPassLFP(LFPdata,Fs,lowCutoff,highCutoff,alignInit)
% Band-pass filters LFPdata between lowCutoff and highCutoff
% 
% LFPdata is LFP data organized by sweeps (sweeps are rows, different
% samples are columns)
% 
% Fs is current sampling frequency of LFPdata
% 
% lowCutoff is lower cut-off for band-pass filter in Hz
% highCutoff is upper cut-off for band-pass filter in Hz
%
% if alignInit is 1, LFP traces for each sweep are aligned at their initial
% values (all initial values set to 0, and other LFP values for sweeps modified
% accordingly); else if alignInit is 0, LFP data is raw

% Low-pass-filter data
disp('LP-filtering');
LFPdata=fftFilter(LFPdata',Fs,highCutoff,1);
LFPdata=LFPdata';
disp('Done LP-filtering');

% High-pass-filter data
disp('HP-filtering');
LFPdata=fftFilter(LFPdata',Fs,lowCutoff,2);
LFPdata=LFPdata';
disp('Done HP-filtering');

% Plot band-pass-filtered data
% Align all LFP traces to initial value = 0
if alignInit==1
    for i=1:size(LFPdata,1)
        initVal=LFPdata(i,1);
        LFPdata(i,:)=LFPdata(i,:)-initVal;
    end
end

bandPassedLFP=LFPdata;
