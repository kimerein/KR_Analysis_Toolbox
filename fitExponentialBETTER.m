function [A,tau,normOfResid,subX,subY,forcezero_tau,fitout,f2]=fitExponentialBETTER(x,y,LEDonset,waitS)
% x and y are row vectors
% LEDonset is the time in seconds when the LED comes on
fitout=[];
f2=[];
% Fit to the following function form:
% y=A*exp(-t/tau)

% Wait for first effects
% waitS=0.003; % Dist. across expts. might indicate time to silence LGN/LP
% waitS=0.01;
% finishS=waitS+0.05; % for anesthetized
% finishS=waitS+0.024;
% Awindow=[-0.05 0];
% Awindow=[-0.05 0];
Awindow=[-0.05 0];
% finishS=waitS+0.041;
% finishS=waitS+0.045;
% finishS=waitS+0.05;
% finishS=waitS+0.033;
finishS=waitS+0.06;
% finishS=waitS+0.2;
% finishS=waitS+0.085;

% finishS=waitS+0.125;
% zeroWindow=[0.01 0.02]; % for dLGN
% finishS=waitS+zeroWindow(2);
% zeroWindow=[0.037 0.05]; % for anesthetized
% zeroWindow=[0.03 0.041]; % for anesthetized
% zeroWindow=[0.035 0.045]; % for anesthetized
zeroWindow=[0.045 0.055]; % for anesthetized
% zeroWindow=[0.03 0.035]; % for anesthetized
% zeroWindow=[0.12 0.14]; % for anesthetized
% zeroWindow=[0.195 0.2]; % for anesthetized
% zeroWindow=[0.07 0.085]; % for anesthetized
% zeroWindow=[0.019 0.024];
% zeroWindow=[0.063 0.1];
% zeroWindow=[0.02 0.06];
% zeroWindow=[0.1 0.125];

% finishS=waitS+0.1;
% zeroWindow=[0.08 0.1];

% Fit to y data from LEDonset+0.003 seconds
% to LEDonset+0.05 seconds
% Get this data
subX=x(x>=LEDonset+waitS & x<=LEDonset+finishS);
subX=subX-min(subX);
subY=y(x>=LEDonset+waitS & x<=LEDonset+finishS);
forcedA=mean(y(x>=LEDonset+Awindow(1) & x<=LEDonset+Awindow(2)));

% Zero curve during zeroWindow
subY=subY-mean(subY(subX>=zeroWindow(1) & subX<=zeroWindow(2)));
% subY=subY-min(subY(subX>=zeroWindow(1) & subX<=zeroWindow(2)));
firstZero=find(subY<=0,1,'first');

% fo=fitoptions('Method','NonlinearLeastSquares','StartPoint',[max(subY) 0.01 max(subY)/10 0.1]);
% ft=fittype('a.*exp(-x./b)+c.*exp(-x./d)','options',fo);
% [fitout,gof,output]=fit(subX,subY,ft);
% figure(); plot(subX,subY,'Color','k'); hold on; plot(fitout);
% newx=0:0.0001:0.2;
% hold on; plot(newx,fitout.a.*exp(-newx./fitout.b),'Color','g');
% hold on; plot(newx,fitout.c.*exp(-newx./fitout.d),'Color','b');

figure(); 
f2 = fit(subX,subY,'exp2');
plot(f2,subX,subY);
disp(f2);

% Get log of subY
log_subY=log(subY(1:firstZero-1)); % linear if y is exponential
log_subX=subX(1:firstZero-1);
logForcedA=log(forcedA);

% Fit a line to log_subY
p=polyfit(log_subX',log_subY',1);
% fzt=mean((log_subY(2:end)-log_subY(1))./log_subX(2:end));
fzt=mean((log_subY(2:end)-logForcedA)./log_subX(2:end));
forcezero_tau=-1/fzt;

% Get exponential fit from linear fit
% p(1) is the slope, p(2) is the y-intercept
tau=-1/p(1);
A=exp(p(2));
resid=subY-A*exp(-subX/tau);
normOfResid=(sum(resid.^2))^(1/2);

% Plot results
figure();
plot(subX,subY,'Color','b');
hold on;
plot(subX,A*exp(-subX/tau),'Color','r');
title('Tau with fit A');

figure();
plot(subX,subY,'Color','b');
hold on;
plot(subX,subY(1)*exp(-subX/forcezero_tau),'Color','r');
title('Tau with force A');
