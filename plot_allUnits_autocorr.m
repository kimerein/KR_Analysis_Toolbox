function plot_allUnits_autocorr(spikes,useAssigns)

useLEDcond=60;

for i=1:length(useAssigns)
    a=useAssigns(i);
    [ev,sp,useLEDcond,maxlag]=plot_unit_autocorr(spikes,a,useLEDcond);
    allEv.x(i,:)=ev.x;
    allEv.y(i,:)=ev.y;
    allSp.x(i,:)=sp.x;
    allSp.y(i,:)=sp.y;
end

figure(); 
m=mean(allEv.y,1);
ymax=max(m);
plot(allEv.x(1,:),m,'Color','k');
hold on;
for i=0:(1/useLEDcond)*1000:maxlag*1000
    if i==0
    else
        line([i i],[0 ymax],'Color','red');
        line([-i -i],[0 ymax],'Color','red');
    end
end
title('Evoked');

figure(); 
m=mean(allSp.y,1);
ymax=max(m);
plot(allSp.x(1,:),m,'Color','k');
hold on; 
for i=0:(1/useLEDcond)*1000:maxlag*1000
    if i==0
    else
        line([i i],[0 ymax],'Color','red');
        line([-i -i],[0 ymax],'Color','red');
    end
end
title('Spontaneous');