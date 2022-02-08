function out=getShiftWithLED(psth,thetaDiff,noTheta_con,noTheta_led,ledOffvals,ledOnvals,stimWindow,spontWindow,useTheseStim)

if isstruct(psth)
    l=psth.unitLED{1};
else
    l=psth;
end
[n,x]=hist(nanmean(thetaDiff(ismember(l,ledOffvals),:),2),-0.2+0.002:0.004:0.2);
figure(); plot(x,n,'Color','k'); hold on;
[n,x]=hist(nanmean(thetaDiff(ismember(l,ledOnvals),:),2),-0.2+0.002:0.004:0.2);
plot(x,n,'Color','c');
conthetadiff=nanmedian(nanmean(thetaDiff(ismember(l,ledOffvals),:),2));
ledthetadiff=nanmedian(nanmean(thetaDiff(ismember(l,ledOnvals),:),2));

if isempty(useTheseStim)
    useStim=1:length(noTheta_con);
else
    useStim=useTheseStim;
end
for i=useStim
    t=noTheta_con(1).allS.t;
    f=noTheta_con(1).allS.f;
    if isempty(spontWindow)
        F1con(i)=nanmean(nanmean(noTheta_con(i).F1amp(:,t>stimWindow(1) & t<=stimWindow(2)),1),2);
        F1led(i)=nanmean(nanmean(noTheta_led(i).F1amp(:,t>stimWindow(1) & t<=stimWindow(2)),1),2);
    else
        F1con(i)=nanmean(nanmean(noTheta_con(i).F1amp(:,t>stimWindow(1) & t<=stimWindow(2)),1),2) - nanmean(nanmean(noTheta_con(i).F1amp(:,t>spontWindow(1) & t<spontWindow(2)),1),2);
        F1led(i)=nanmean(nanmean(noTheta_led(i).F1amp(:,t>stimWindow(1) & t<=stimWindow(2)),1),2) - nanmean(nanmean(noTheta_led(i).F1amp(:,t>spontWindow(1) & t<spontWindow(2)),1),2);
    end
    if i==1
        sumF1con=nanmean(noTheta_con(i).F1amp,1);
        sumF1led=nanmean(noTheta_led(i).F1amp,1);
    else
        sumF1con=sumF1con+nanmean(noTheta_con(i).F1amp,1);
        sumF1led=sumF1led+nanmean(noTheta_led(i).F1amp,1);
    end
end

figure(); 
plot(t,sumF1con,'Color','k'); 
hold on;
plot(t,sumF1led,'Color','c'); 
line([stimWindow(1) stimWindow(1)],[nanmin(sumF1con) nanmax(sumF1con)],'Color','r');
line([stimWindow(2) stimWindow(2)],[nanmin(sumF1con) nanmax(sumF1con)],'Color','r');

disp([conthetadiff; ledthetadiff; conthetadiff-ledthetadiff; nanmean(F1led./F1con)]);
out=nanmean(F1led./F1con);

end