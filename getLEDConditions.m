function [ledConds,mForLED]=getLEDConditions(LEDbySweep,Fs)
% This function may need some upgrading for new LED configuration on rig
% For very low LED voltages, may be problems with this function
% Works to pick out LED trials from the recorded LED data, rather than
% trusting the saved LED conditions file, which sometimes says the wrong
% trials had LED on
%
% RETURNS:
% ledConds: a vector with 0's for trials with LED off and 1's for trials
% with LED on
% mForLED: the maximum LED signal
%
% PARAMETERS:
% LEDbySweep: a matrix with samples as columns and trials as rows that
% contains data acquired on the LED channel (feed LED output directly back
% into daq board)
% Fs: sampling rate of LEDbySweep

global saveToDir2

ledConds=zeros(1,size(LEDbySweep,1));
thisFig=figure;
title('Check LED detection');
for i=1:size(LEDbySweep,1)
    m=max(LEDbySweep(i,:));
    if m>0.003
        ledConds(i)=1;
    else
        f=fftFilter(LEDbySweep(i,:)',Fs,15,1);
        if max(f)>0.0003
            ledConds(i)=1;
        else
            ledConds(i)=0;
        end
        m=max(f);
    end
    hold on;
    if ledConds(i)==1
        plot(LEDbySweep(i,:),'Color','r');
    else
        plot(LEDbySweep(i,:),'Color','k');
    end
    if i==1 
        mForLED=m;
    else
        if m>mForLED
            mForLED=m;
        end
    end
end
saveas(thisFig,strcat(saveToDir2,'\checkLEDdetect.fig'));
