% Script for checking LED saved conditions

figure();
for i=1:length(allLEDConds)
    if allLEDConds(i)==0
        plot(ledData(i,:),'Color','black');
    else
        plot(ledData(i,:),'Color','red');
    end
    hold on;
end