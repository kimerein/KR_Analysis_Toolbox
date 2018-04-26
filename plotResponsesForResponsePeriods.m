

responseWindows=[0.5 2;2 2.05; 2.05 2.08; 2.08 2.1; 2.1 2.2; 2.2 2.4; 2.4 2.6; 2.6 2.8; 2.8 3; 3 3.2];

% Seconds from 2.06
timesFromResponseOnset=[0;0.01;0.04;0.08;0.12;0.24;0.48;1.92];
timesFromResponseOnset=timesFromResponseOnset-0.04;

responses=[-0.0102875;-0.0228125;-0.0428625;-0.020575;0.0072875;-0.0001375;-0.00255;-0.01585];

figure;
title('Averages of Whole Response Period');
plot(timesFromResponseOnset,responses,'-or');

partResponses=[-0.031875;-0.040625;-0.067525;-0.0283;0.00325;-0.00035;0.00845;-0.0138];
figure;
title('Averages of Cols. D through G');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.024;-0.03128;-0.05656;-0.0258;0.00994;-0.00224;0.00574;-0.01414];
figure;
title('Averages of Cols. D through H');
plot(timesFromResponseOnset,partResponses,'-or');

part2Responses=[-0.0441;-0.0559;-0.0427;-0.011;0.019466667;0.000833333;-0.0039;-0.0102];
figure;
title('Averages of Only After Initial Dip - Pink');
plot(timesFromResponseOnset,part2Responses,'-or');

partResponses=[-0.060575;-0.0597;-0.0576;-0.0283;0.00325;-0.00035;0.007;0.0155];
figure;
title('Averages of Cols. D through G for all Aligned to Baseline');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.05276;-0.05036;-0.04664;-0.0258;0.00994;-0.00224;0.0043;0.01516];
figure;
title('Averages of Cols. D through H for all Aligned to Baseline');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.072933333;-0.074966667;-0.032766667;-0.011;0.019466667;0.000833333;-0.005333333;0.0192];
figure;
title('Averages of Cols. F through H for all Aligned to Baseline');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.0634;-0.06344;-0.02834;-0.01158;0.02066;-0.00014;-0.00698;0.01794];
figure;
title('Averages of Cols. F through H for all Aligned to Baseline -WEIGHTED');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.00973;-0.0274;-0.03116;-0.01291;0.01205;0.00195;-0.0126;-0.01696];
figure;
title('Averages of Cols. F through all WEIGHTED');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.017032;-0.022202;-0.021432;-0.006412;0.010842;9.2E-05;-0.002304;-0.006448];
figure;
title('Averages of Cols. E through H for all WEIGHTED');
plot(timesFromResponseOnset,partResponses,'-or');

partResponses=[-0.017032;-0.022202;-0.021432;-0.006412;0.010842;9.2E-05;-0.002304;0.00884];
figure;
title('Averages of Cols. E through H for all WEIGHTED 2');
plot(timesFromResponseOnset,partResponses,'-or');


% Curve fitting 
figure;
partResponses=[-0.022202;-0.021432;-0.006412;0.010842;9.2E-05;-0.002304;-0.006448];
timesFromResponseOnset2=timesFromResponseOnset(2:end);
%timesFromResponseOnset2=timesFromResponseOnset2-timesFromResponseOnset2(1);
partResponses2=partResponses-min(partResponses);
partResponses2=max(partResponses2)-partResponses2;
fit=polyfit(timesFromResponseOnset2,log(partResponses2),1);
%semilogy(timesFromResponseOnset(2:end),partResponses,'o',timesFromResponseOnset(2:end),exp(fit(2)).*exp(fit(1)*timesFromResponseOnset2));
semilogy(timesFromResponseOnset(2:end),partResponses2,'o',timesFromResponseOnset(2:end),exp(fit(2)).*exp(fit(1)*timesFromResponseOnset2));

figure;
r=[0.010842;9.2E-05;-0.002304;-0.006448];
r=r-min(r)+0.0001;
plot(timesFromResponseOnset(end-3:end),r);

figure;
plot(timesFromResponseOnset(end-3:end),log(r));