

x=0:0.0012:3;
ystart=chirp(x,1,3,30,'logarithmic');
y=[ystart fliplr(ystart) ystart];
figure(); 
subplot(3,1,1);
neuralx=linspace(0,10.5,length(downSampAv(nanmean(n_V1,1),20)));
plot(neuralx,downSampAv(nanmean(n_V1,1),20));
base=nanmean(nanmean(n_V1(:,neuralx<=0.5),1),2);
hold on; plot([0 10.5],[base base],'Color','g');
xlim([0 10.5]);
subplot(3,1,2);
plot(linspace(0,9,length(y))+0.5,y,'Color','k');
xlim([0 10.5]);
subplot(3,1,3);
xforchirp=linspace(0,3,length(ystart));
chirpfreq=1*((30/1)^(1./3)).^xforchirp;
plot(0.5+linspace(0,9,length(y)),[chirpfreq fliplr(chirpfreq) chirpfreq]);
xlim([0 10.5]);