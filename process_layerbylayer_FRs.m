n=47;

for i=1:16
curr_noLED=layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{1};
noLED_means_spont(i)=mean(curr_noLED(1:3:n));
LED_means_spont(i)=mean(curr_LED(1:3:n));
noLED_stds_spont(i)=std(curr_noLED(1:3:n));
LED_stds_spont(i)=std(curr_LED(1:3:n));
end
for i=1:16
curr_noLED=layers{i}.FRs_noLED{2}-layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{2}-layers{i}.FRs_LED{1};
noLED_means_ev(i)=mean(curr_noLED(1:3:n));
LED_means_ev(i)=mean(curr_LED(1:3:n));
noLED_stds_ev(i)=std(curr_noLED(1:3:n));
LED_stds_ev(i)=std(curr_LED(1:3:n));
end
noLED_means_spont_lines{1}=noLED_means_spont;
LED_means_spont_lines{1}=LED_means_spont;
noLED_stds_spont_lines{1}=noLED_stds_spont;
LED_stds_spont_lines{1}=LED_stds_spont;
noLED_means_ev_lines{1}=noLED_means_ev;
LED_means_ev_lines{1}=LED_means_ev;
noLED_stds_ev_lines{1}=noLED_stds_ev;
LED_stds_ev_lines{1}=LED_stds_ev;

for i=1:16
curr_noLED=layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{1};
noLED_means_spont(i)=mean(curr_noLED(2:3:n));
LED_means_spont(i)=mean(curr_LED(2:3:n));
noLED_stds_spont(i)=std(curr_noLED(2:3:n));
LED_stds_spont(i)=std(curr_LED(2:3:n));
end
for i=1:16
curr_noLED=layers{i}.FRs_noLED{2}-layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{2}-layers{i}.FRs_LED{1};
noLED_means_ev(i)=mean(curr_noLED(2:3:n));
LED_means_ev(i)=mean(curr_LED(2:3:n));
noLED_stds_ev(i)=std(curr_noLED(2:3:n));
LED_stds_ev(i)=std(curr_LED(2:3:n));
end
noLED_means_spont_lines{2}=noLED_means_spont;
LED_means_spont_lines{2}=LED_means_spont;
noLED_stds_spont_lines{2}=noLED_stds_spont;
LED_stds_spont_lines{2}=LED_stds_spont;
noLED_means_ev_lines{2}=noLED_means_ev;
LED_means_ev_lines{2}=LED_means_ev;
noLED_stds_ev_lines{2}=noLED_stds_ev;
LED_stds_ev_lines{2}=LED_stds_ev;

for i=1:16
curr_noLED=layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{1};
noLED_means_spont(i)=mean(curr_noLED(3:3:n));
LED_means_spont(i)=mean(curr_LED(3:3:n));
noLED_stds_spont(i)=std(curr_noLED(3:3:n));
LED_stds_spont(i)=std(curr_LED(3:3:n));
end
for i=1:16
curr_noLED=layers{i}.FRs_noLED{2}-layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{2}-layers{i}.FRs_LED{1};
noLED_means_ev(i)=mean(curr_noLED(3:3:n));
LED_means_ev(i)=mean(curr_LED(3:3:n));
noLED_stds_ev(i)=std(curr_noLED(3:3:n));
LED_stds_ev(i)=std(curr_LED(3:3:n));
end
noLED_means_spont_lines{3}=noLED_means_spont;
LED_means_spont_lines{3}=LED_means_spont;
noLED_stds_spont_lines{3}=noLED_stds_spont;
LED_stds_spont_lines{3}=LED_stds_spont;
noLED_means_ev_lines{3}=noLED_means_ev;
LED_means_ev_lines{3}=LED_means_ev;
noLED_stds_ev_lines{3}=noLED_stds_ev;
LED_stds_ev_lines{3}=LED_stds_ev;

 figure(); plot(1:16,noLED_means_spont_lines{1},'Color','k'); hold on;
 plot(1:16,noLED_means_spont_lines{2},'Color','k'); hold on;
 plot(1:16,noLED_means_spont_lines{3},'Color','k');
hold on;
plot(1:16,LED_means_spont_lines{1},'Color','r'); hold on;
 plot(1:16,LED_means_spont_lines{2},'Color','r'); hold on;
 plot(1:16,LED_means_spont_lines{3},'Color','r');
 
 figure(); plot(1:16,noLED_means_ev_lines{1},'Color','k'); hold on;
 plot(1:16,noLED_means_ev_lines{2},'Color','k'); hold on;
 plot(1:16,noLED_means_ev_lines{3},'Color','k');
hold on;
plot(1:16,LED_means_ev_lines{1},'Color','r'); hold on;
 plot(1:16,LED_means_ev_lines{2},'Color','r'); hold on;
 plot(1:16,LED_means_ev_lines{3},'Color','r');
 
 for i=1:16
curr_noLED=layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{1};
noLED_means_spont(i)=mean(curr_noLED);
LED_means_spont(i)=mean(curr_LED);
noLED_stds_spont(i)=std(curr_noLED);
LED_stds_spont(i)=std(curr_LED);
end
for i=1:16
curr_noLED=layers{i}.FRs_noLED{2}-layers{i}.FRs_noLED{1};
curr_LED=layers{i}.FRs_LED{2}-layers{i}.FRs_LED{1};
noLED_means_ev(i)=mean(curr_noLED);
LED_means_ev(i)=mean(curr_LED);
noLED_stds_ev(i)=std(curr_noLED);
LED_stds_ev(i)=std(curr_LED);
end

 figure(); errorbar(1:16,noLED_means_spont,noLED_stds_spont,'Color','k'); hold on; errorbar(1:16,LED_means_spont,LED_stds_spont,'Color','r');
 figure(); errorbar(1:16,noLED_means_ev,noLED_stds_ev,'Color','k'); hold on; errorbar(1:16,LED_means_ev,LED_stds_ev,'Color','r');