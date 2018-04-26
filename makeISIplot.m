function [diffgamma,diffgammaNorm,x,y]=makeISIplot(spikes,binsize)

% timeBins={[0 4]; [6.5 10]};
timeBins={[0 2.5]; [8 10]};
% timeBins={[0 2.5]};
% timeBins={[0 10]};
% timeBins={[0 2]};

spiketimes=spikes.unwrapped_times;
takeSpiketimes=[];
for i=1:length(timeBins)
    tb=timeBins{i};
    disp(tb);
    ts=spiketimes(spikes.spiketimes>=tb(1) & spikes.spiketimes<=tb(2));
    takeSpiketimes=[takeSpiketimes ts];
end

spikes_forfig.unwrapped_times=sort(takeSpiketimes);

figure();
h=axes();
% [~,~,y,x]=plotAutoCorr(spikes_forfig,h,100,binsize);
[x,y]=plot_isi(spikes_forfig.unwrapped_times);
figure(); 
plot([-x x],[smooth(y,3); smooth(y,3)]);

% xout=0:0.001:0.1;
% [n]=histc(diff(takeSpiketimes),xout);
% figure();
% x=[sort(-xout) xout(2:end)];
% y=[fliplr(n) n(2:end)];
% % plot(xout,n);
% plot(x,y);
% title('Histogram of ISIs');
% figure();
% h=axes();
% % bar(h,xout,n);
% bar(h,x,y);
% % x=xout; y=n;

params.Fs=(1/(x(2)-x(1)))*1000;
params.tapers=[3 10];
params.fpass=[30 80];
params.trialave=1;
[S,f]=mtspectrumpb(y',params);

% figure();
% plot(f,S);
% title('frequency content of ISI');

% Is there a narrowband peak?
narrowgamma=nanmean(S(f>=55 & f<=65));
broadgamma=nanmean(S(f<=55 | f>=65));
diffgamma=narrowgamma-broadgamma;

% Is there a narrowband peak?
S=S-min(S);
S=S./max(S);
narrowgamma=nanmean(S(f>=55 & f<=65));
broadgamma=nanmean(S(f<=55 | f>=65));
diffgammaNorm=narrowgamma-broadgamma;

end

function [x,n]=plot_isi(spiketimes)
       
    data.isi_maxlag=0.1;

    maxlag = data.isi_maxlag;
    bins=100;
    
    % make plot
    isis = diff(spiketimes);
    isis = isis(isis <= maxlag);
    [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
    ymax   = max(n)+1;
    
    b2 = bar(x,n,1.0); hold off
    set(b2,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])
    hold on
    b2 = bar(-x,n,1.0); hold off
    set(b2,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])
    
    % update axes
    set(gca,'YLim',[0 ymax],'XLim',[-1000*maxlag 1000*maxlag])
    xlabel('Interspike interval (msec)')
    ylabel('No. of spikes')

end
