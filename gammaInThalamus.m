function out=gammaInThalamus(spikes,useAssigns)

% ledOn=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60]+0.05;
% ledOff=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% ledOn=[1 2 4 6 8 10]+0.05;
% ledOff=[1 2 4 6 8 10];
% ledOn=0.05;
% ledOff=0;
ledOn=5;
ledOff=0;

[~,~,~,~,~,~,~,~,~,~,~,~,forfigs]=plot_unit_autocorr(spikes,useAssigns,ledOff);
[~,~,~,~,~,~,~,~,~,~,~,~,forfigs_led]=plot_unit_autocorr(spikes,useAssigns,ledOn);

figure(); 
plot(forfigs.evoked_autocorr.x,forfigs.evoked_autocorr.y,'Color','k');
hold on; 
plot(forfigs_led.evoked_autocorr.x,forfigs_led.evoked_autocorr.y,'Color','r');
title('Auto-Corr With and Without LED');
% xlim([-100 100]);
xlim([-40 40]);
ylim([0 max([forfigs.evoked_autocorr.y forfigs_led.evoked_autocorr.y])]);
% mi=min([forfigs.evoked_autocorr.y(forfigs.evoked_autocorr.x<=-10 & forfigs.evoked_autocorr.x>=-100) forfigs_led.evoked_autocorr.y(forfigs_led.evoked_autocorr.x<=-10 & forfigs_led.evoked_autocorr.x>=-100)]);
% ma=max([forfigs.evoked_autocorr.y forfigs_led.evoked_autocorr.y]);
% ylim([mi ma]);
out.autocorr_x=forfigs.evoked_autocorr.x;
out.autocorr_y=forfigs.evoked_autocorr.y;
out.autocorr_x_led=forfigs_led.evoked_autocorr.x;
out.autocorr_y_led=forfigs_led.evoked_autocorr.y;

figure(); 
plot(forfigs.evoked_autocorr_power.x,forfigs.evoked_autocorr_power.y,'Color','k');
hold on; 
plot(forfigs_led.evoked_autocorr_power.x,forfigs_led.evoked_autocorr_power.y,'Color','r');
title('Power Spec of Auto-Corr');
xlim([1 250]);
ma=max([forfigs.evoked_autocorr_power.y(forfigs.evoked_autocorr_power.x>=4); forfigs_led.evoked_autocorr_power.y(forfigs.evoked_autocorr_power.x>=4)]);
ylim([0 ma]);

figure(); 
plot(forfigs.evoked_autocorr_power.x,forfigs.evoked_autocorr_power.y,'Color','k');
hold on; 
plot(forfigs_led.evoked_autocorr_power.x,forfigs_led.evoked_autocorr_power.y,'Color','r');
title('Power Spec of Auto-Corr');
xlim([30 80]);
out.power_x=forfigs.evoked_autocorr_power.x;
out.power_y=forfigs.evoked_autocorr_power.y;
out.power_x_led=forfigs_led.evoked_autocorr_power.x;
out.power_y_led=forfigs_led.evoked_autocorr_power.y;

