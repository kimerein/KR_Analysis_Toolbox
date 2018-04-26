function [A,tau,normOfResid,subX,subY,forcezero_tau]=fitExponential(x,y,LEDonset,waitS)
% x and y are row vectors
% LEDonset is the time in seconds when the LED comes on

% Fit to the following function form:
% y=A*exp(-t/tau)

% Wait for first effects
% waitS=0.003; % Dist. across expts. might indicate time to silence LGN/LP
% waitS=0.01;
% finishS=waitS+0.05; % for anesthetized
% finishS=waitS+0.024;
% Awindow=[-0.05 0];
% Awindow=[-0.05 0];
Awindow=[-0.1 0];
finishS=waitS+0.033;

% finishS=waitS+0.125;
% zeroWindow=[0.01 0.02]; % for dLGN
% finishS=waitS+zeroWindow(2);
% zeroWindow=[0.037 0.05]; % for anesthetized
zeroWindow=[0.03 0.033];
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
firstZero=find(subY<=0,1,'first');

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
