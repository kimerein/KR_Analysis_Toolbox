downTo=0; % where 1000 is 100% of response remaining

tau_x=x;
tau_y=avled;

bottomVal=mean(tau_y(tau_x>50 & tau_x<100));
tau_y=tau_y-bottomVal;
topVal=mean(tau_y(tau_x>-20 & tau_x<0));
tau_y=tau_y*((1000-downTo)/topVal);
tau_y=tau_y+downTo;

plot(downSampAv(tau_x,2),downSampAv(tau_y,2));
hold all;