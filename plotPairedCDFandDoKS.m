function plotPairedCDFandDoKS(nameOfFiles,t)

a=load([nameOfFiles '_noThetaNoLED.mat']);
temp=a.temp;
evoked=nanmean(temp(:,t>4 & t<6.5),2)-nanmean(temp(:,t>3 & t<3.75),2);
noThetaNoLED_evoked=evoked;

a=load([nameOfFiles '_thetaNoLED.mat']);
temp=a.temp;
evoked=nanmean(temp(:,t>4 & t<6.5),2)-nanmean(temp(:,t>3 & t<3.75),2);
thetaNoLED_evoked=evoked;

a=load([nameOfFiles '_noThetaLED.mat']);
temp=a.temp;
evoked=nanmean(temp(:,t>4 & t<6.5),2)-nanmean(temp(:,t>3 & t<3.75),2);
noThetaLED_evoked=evoked;

a=load([nameOfFiles '_thetaLED.mat']);
temp=a.temp;
evoked=nanmean(temp(:,t>4 & t<6.5),2)-nanmean(temp(:,t>3 & t<3.75),2);
thetaLED_evoked=evoked;

currdiff=thetaNoLED_evoked-noThetaNoLED_evoked; 
figure(); 
[f,x,flo,fup]=ecdf(currdiff,'Bounds','on');
plot(x,f,'Color','k'); hold on;
plot(x,flo,'Color','k'); 
plot(x,fup,'Color','k');
test_cdf=makedist('tlocationscale','mu',nanmean(currdiff),'sigma',std(currdiff,[],1),'nu',1);
y=cdf(test_cdf,x);
plot(x,y,'Color','r');
%x_values=linspace(min(currdiff),max(currdiff));
%plot(x_values,normcdf(x_values,nanmean(currdiff),std(currdiff,[],1)),'r-');
[~,p]=kstest(currdiff,'CDF',test_cdf);
%[~,p]=kstest((currdiff-nanmean(currdiff))./std(currdiff,[],1));
disp('theta no LED vs no theta no LED, ks test p-val: ');
disp(p);

currdiff=thetaLED_evoked-thetaNoLED_evoked; 
figure(); 
ecdf(currdiff,'Bounds','on');
hold on;
x_values=linspace(min(currdiff),max(currdiff));
plot(x_values,normcdf(x_values,nanmean(currdiff),std(currdiff,[],1)),'r-');
[~,p]=kstest((currdiff-nanmean(currdiff))./std(currdiff,[],1));
disp('theta no LED vs theta LED, ks test p-val: ');
disp(p);


end