function [noLED_means_spont,LED_means_spont]=plotLayerDependencies(layers,block)

for i=1:16
    curr_noLED=layers{i}.FRs_noLED{block};
    curr_LED=layers{i}.FRs_LED{block};
    noLED_means_spont(i)=mean(curr_noLED);
    LED_means_spont(i)=mean(curr_LED);
    noLED_stds_spont(i)=std(curr_noLED);
    LED_stds_spont(i)=std(curr_LED);
end
figure(); 
plot(1:16,noLED_means_spont,'Color','k');
hold on;
plot(1:16,LED_means_spont,'Color','r');