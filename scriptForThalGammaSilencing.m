
ledOff=0;
ledOn=5.05;

[~,~,~,~,~,n_V1]=psth_wStd_trialByTrial(criterionSpikes,0.0012*1000,ledOff,ledOn,length(unique(criterionSpikes.trials)),unique(criterionSpikes.trials));

[~,trialInd]=unique(criterionSpikes.trials);
l=criterionSpikes.led(trialInd);

% SMALL WINDOW
% paramsCoh.tapers=[0.05 2];
% [SLED,fLED]=mtspectrumpb(n_V1(floor(l)==floor(ledOn),xtimes>=1.17 & xtimes<=1.375)',paramsCoh);
% figure(); plot(fLED(fLED>=20 & fLED<=100),SLED(fLED>=20 & fLED<=100),'Color','r');
% [S,f]=mtspectrumpb(n_V1(floor(l)==floor(ledOff),xtimes>=1.17 & xtimes<=1.375)',paramsCoh);
% hold on; plot(f(f>=20 & f<=100),S(f>=20 & f<=100),'Color','k');
% BIG WINDOW
paramsCoh.tapers=[5 9];
[SLED,fLED]=mtspectrumpb(n_V1(floor(l)==floor(ledOn),xtimes>=1.2 & xtimes<=2.2)',paramsCoh);
figure(); plot(fLED(fLED>=20 & fLED<=100),SLED(fLED>=20 & fLED<=100),'Color','r');
[S,f]=mtspectrumpb(n_V1(floor(l)==floor(ledOn),xtimes>=2.2 & xtimes<=3.2)',paramsCoh);
hold on; plot(f(f>=20 & f<=100),S(f>=20 & f<=100),'Color','k');
% curr=n_V1(floor(l)==floor(ledOn),xtimes>=1.2 & xtimes<=1.35);
% clear veccurr
% j=1;
% for i=1:size(curr,1)
%     if i+9>size(curr,1)
%         break
%     else
%         c=curr(i:i+9,:);
%     end
%     veccurr(j,:)=c(1:end);
%     j=j+1;
% end
% [SLED,fLED]=mtspectrumpb(veccurr',paramsCoh);
% figure(); plot(fLED(fLED>=20 & fLED<=100),SLED(fLED>=20 & fLED<=100),'Color','r');
% curr=n_V1(floor(l)==floor(ledOff),xtimes>=1.2 & xtimes<=1.35);
% clear veccurr
% j=1;
% for i=1:size(curr,1)
%     if i+9>size(curr,1)
%         break
%     else
%         c=curr(i:i+9,:);
%     end
%     veccurr(j,:)=c(1:end);
%     j=j+1;
% end
% [S,f]=mtspectrumpb(veccurr',paramsCoh);
% hold on; plot(f(f>=20 & f<=100),S(f>=20 & f<=100),'Color','k');

% SMALL WINDOW
% paramsCoh.tapers=[5 9];
% [C5,t5,f5]=mtspecgrampb(n_V1(floor(l)==floor(ledOn),:)',[0.5 0.01],paramsCoh);
% [C0,t0,f0]=mtspecgrampb(n_V1(floor(l)==floor(ledOff),:)',[0.5 0.01],paramsCoh);
% figure(); imagesc(t0,f0(f0>=20 & f0<=100),C0(:,f0>=20 & f0<=100)');
% figure(); imagesc(t5,f5(f5>=20 & f5<=100),C5(:,f5>=20 & f5<=100)');
% BIG WINDOW
paramsCoh.tapers=[5 9];
[C5,t5,f5]=mtspecgrampb(n_V1(floor(l)==floor(ledOn),:)',[1 0.01],paramsCoh);
[C0,t0,f0]=mtspecgrampb(n_V1(floor(l)==floor(ledOff),:)',[1 0.01],paramsCoh);
figure(); imagesc(t0,f0(f0>=20 & f0<=100),C0(:,f0>=20 & f0<=100)');
figure(); imagesc(t5,f5(f5>=20 & f5<=100),C5(:,f5>=20 & f5<=100)');

out.SLED=SLED;
out.fLED=fLED;
out.S=S;
out.f=f;
out.C5=C5;
out.t5=t5;
out.f5=f5;
out.C0=C0;
out.t0=t0;
out.f0=f0;

save('W:\Analysis Computer\New Thalamic Gamma 150429\Silencing Thalamus\mawake96_bigwindow.mat','out');