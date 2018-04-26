%% Baseline
% params.ledCondNum=[4 8];
% params.ledCondDen=[3 5];
params.ledCondNum=[4.04 8.08];
params.ledCondDen=[3.03 5.05];
% params.ledCondNum=[2 6];
% params.ledCondDen=[1 7];
% params.ledCondNum=[2.02 6.06];
% params.ledCondDen=[1.01 7.07];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[0.743 0.99];
params.window=[0.783 1.03];

%% Persistent window
% params.ledCondNum=[4 8];
% params.ledCondDen=[3 5];
params.ledCondNum=[4.04 8.08];
params.ledCondDen=[3.03 5.05];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[1.303 1.55];
params.window=[1.343 1.59];

%% Persistent window but no thalamic silencing
% params.ledCondNum=[2 6];
% params.ledCondDen=[1 7];
params.ledCondNum=[2.02 6.06];
params.ledCondDen=[1.01 7.07];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[1.303 1.55];
params.window=[1.343 1.59];

%% Preceding window
params.ledCondNum=[2 6];
params.ledCondDen=[1 7];
% params.ledCondNum=[2.02 6.06];
% params.ledCondDen=[1.01 7.07];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[1 4];
params.window=[1.092 1.339];

%% Evoked response
% params.ledCondNum=[2 6];
% params.ledCondDen=[1 7];
% params.ledCondNum=[2.02 6.06];
% params.ledCondDen=[1.01 7.07];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[1 4];
params.window=[1.04 4.04];

%% Disinhibition -- spontaneous
% params.ledCondNum=[2 4 6 8];
% params.ledCondDen=[1 3 5 7];
params.ledCondNum=[2.02 4.04 6.06 8.08];
params.ledCondDen=[1.01 3.03 5.05 7.07];
params.stimCondNum=1:32;
params.stimCondDen=1:32;
% params.window=[0.5 1];
params.window=[0.54 1.04];

%% All trode analysis -- persistent no LED
% all_led26.persistentNoLED=s_ledNum;
% all_led17.persistentNoLED=s_ledDen;

all_led26.persistentNoLED=[all_led26.persistentNoLED s_ledNum];
all_led17.persistentNoLED=[all_led17.persistentNoLED s_ledDen];

%% All trode analysis -- other baseline
% all_led26.baseline=s_ledNum;
% all_led17.baseline=s_ledDen;

all_led26.baseline=[all_led26.baseline s_ledNum];
all_led17.baseline=[all_led17.baseline s_ledDen];

%% All trode analysis --  baseline
% all_led48.baseline=s_ledNum;
% all_led35.baseline=s_ledDen;

all_led48.baseline=[all_led48.baseline s_ledNum];
all_led35.baseline=[all_led35.baseline s_ledDen];

%% All trode analysis --  persistentWindow
% all_led48.persistentWindow=s_ledNum;
% all_led35.persistentWindow=s_ledDen;

all_led48.persistentWindow=[all_led48.persistentWindow s_ledNum];
all_led35.persistentWindow=[all_led35.persistentWindow s_ledDen];
%% All trode analysis
all_led1357.disinhibition=led1357.disinhibition;
all_led17.evokedResponse=led17.evokedResponse;
all_led2468.disinhibition=led2468.disinhibition;
all_led26.evokedResponse=led26.evokedResponse;
all_led35.baseline=led35.baseline;
all_led35.persistentWindow=led35.persistentWindow;
all_led48.baseline=led48.baseline;
all_led48.persistentWindow=led48.persistentWindow;

% all_led1357.disinhibition=[all_led1357.disinhibition led1357.disinhibition];
% all_led17.evokedResponse=[all_led17.evokedResponse led17.evokedResponse];
% all_led2468.disinhibition=[all_led2468.disinhibition led2468.disinhibition];
% all_led26.evokedResponse=[all_led26.evokedResponse led26.evokedResponse];
% all_led35.baseline=[all_led35.baseline led35.baseline];
% all_led35.persistentWindow=[all_led35.persistentWindow led35.persistentWindow];
% all_led48.baseline=[all_led48.baseline led48.baseline];
% all_led48.persistentWindow=[all_led48.persistentWindow led48.persistentWindow];

%% Figs
% figure(); scatter(all_led2468.disinhibition-all_led1357.disinhibition,(all_led48.persistentWindow-all_led48.baseline)./(all_led26.evokedResponse-all_led48.baseline));
% title('spont disinhibition vs persistent as fraction of evoked');
% figure(); scatter(all_led26.evokedResponse-all_led17.evokedResponse,(all_led48.persistentWindow-all_led48.baseline)./(all_led26.evokedResponse-all_led48.baseline));
% title('ev disinhibition vs persistent as fraction of evoked');
% figure(); scatter(all_led17.evokedResponse,(all_led48.persistentWindow-all_led48.baseline)./(all_led26.evokedResponse-all_led48.baseline));
% title('ev response vs persistent as fraction of evoked');
% figure(); scatter(all_led2468.disinhibition-all_led1357.disinhibition,(all_led48.persistentWindow-all_led48.baseline)-(all_led35.persistentWindow-all_led35.baseline));
% title('spont disinhibition vs persistent as amber magnitude minus control');
% figure(); scatter(all_led26.evokedResponse-all_led17.evokedResponse,(all_led48.persistentWindow-all_led48.baseline)./(all_led26.evokedResponse-all_led48.baseline));
% title('ev disinhibition vs persistent as amber magnitude minus control');
% figure(); scatter(all_led17.evokedResponse,(all_led48.persistentWindow-all_led48.baseline)./(all_led26.evokedResponse-all_led48.baseline));
% title('ev response vs persistent as amber magnitude minus control');


figure(); 
% a=all_led17.evokedResponse-all_led35.baseline;
% b=all_led26.evokedResponse-all_led48.baseline;
% a=all_led17.precedingWindow-all_led35.baseline;
% b=all_led26.precedingWindow-all_led48.baseline;
a=all_led17.persistentNoLED-all_led17.baseline;
b=all_led26.persistentNoLED-all_led26.baseline;
c=all_led35.persistentWindow-all_led35.baseline;
d=all_led48.persistentWindow-all_led48.baseline;
% d=c.*(b./a);
useUnits=(a>0) & (b>0);
disp('total # units');
disp(length(a));
disp('Number of units with a and b not equal to zero');
disp(sum(useUnits));

figure(); scatter(c(useUnits)./a(useUnits),d(useUnits)./b(useUnits));

%% 
% p=(d(useUnits)./b(useUnits))-(c(useUnits)./a(useUnits));
% p=(d(useUnits)./b(useUnits))./(c(useUnits)./a(useUnits));
useUnits=(a>0) & (b>0) & (b>d) & (a>c) & (c>0) & (d>0);
disp('total # units');
disp(length(a));
disp('Number of units with all criteria satisfied');
disp(sum(useUnits));
p=(d(useUnits)./b(useUnits))./((d(useUnits)./b(useUnits))+(c(useUnits)./a(useUnits)));
[n,xout]=hist(p(~isinf(p)),20);
plot(xout,n);
title('Hist of Persistence Index');

% figure(); scatter(all_led2468.disinhibition(useUnits)-all_led1357.disinhibition(useUnits),p);
figure(); scatter(all_led2468.disinhibition(useUnits)./all_led1357.disinhibition(useUnits),p);
% figure(); scatter(all_led2468.disinhibition(useUnits)./(all_led2468.disinhibition(useUnits)+all_led1357.disinhibition(useUnits)),p);
title('spont disinhibition vs persistence index');
% figure(); scatter(all_led26.evokedResponse(useUnits)-all_led17.evokedResponse(useUnits),p);
figure(); scatter(b(useUnits)./a(useUnits),p);
title('ev disinhibition vs persistence index');
% figure(); scatter(all_led17.evokedResponse(useUnits),p);
figure(); scatter(all_led17.evokedResponse(useUnits)-all_led35.baseline(useUnits),p);
title('ev response vs persistence index');
figure(); scatter(all_led26.evokedResponse(useUnits)-all_led48.baseline(useUnits),p);
title('amber ev response vs persistence index');