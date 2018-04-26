function [xout1,n1,xout2,n2,FRs1,FRs2]=exploring_layer_dependencies(depend_spikes,ind_spikes)

LEDconds=[3];
LEDwindow=[3.6 4.1];

% Histogram of spiking during LED window for depend spikes
[m s FRs1]=calcMeanAndStdDuringWindow(filtspikes(depend_spikes,0,'led',LEDconds),LEDwindow);
[n1 xout1]=hist(FRs1,10);
% figure(); bar(xout,n);
disp(m); disp(s);

% Histogram of spiking during LED window for ind spikes
[m s FRs2]=calcMeanAndStdDuringWindow(filtspikes(ind_spikes,0,'led',LEDconds),LEDwindow);
[n2 xout2]=hist(FRs2,10);
% figure(); bar(xout,n);
disp(m); disp(s);

figure(); 
subplot(2,1,1); bar(xout1,n1);
subplot(2,1,2); bar(xout2,n2);