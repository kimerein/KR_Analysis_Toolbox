function [f1,S1,f2,S2]=chronuxPowerSpectrum(x,y1,y2,useWindowForSpectrum)

 params.Fs=1/(x(2)-x(1));
%  params.tapers=[0.9 x(end)-x(1) 0];
 params.tapers=[0.9 useWindowForSpectrum(2)-useWindowForSpectrum(1) 0];
 [S1,f1]=mtspectrumpb(y1(:,x>=useWindowForSpectrum(1) & x<=useWindowForSpectrum(2))',params);
 [S2,f2]=mtspectrumpb(y2(:,x>=useWindowForSpectrum(1) & x<=useWindowForSpectrum(2))',params);
 
%  figure(); 
 S1=nanmean(S1,2);
%  plot(f1,S1,'Color','k');
%  axis([0 100 min(S1) max(S1)]);
%  hold on;
 S2=nanmean(S2,2);
%  plot(f2,S2,'Color','r');
%  axis([0 100 min([min(S1) min(S2)]) max([max(S1) max(S2)])]);
 
 return
 
 params.tapers=[0.5 x(end)-x(1) 0];
 plotSpectrogram(y1,[0.5 0.1],params);
 plotSpectrogram(y2,[0.5 0.1],params);
 
end

function plotSpectrogram(y,movingwindow,params)

[S,t,f,R]=mtspecgrampb(y,movingwindow,params);
 subS=S(:,f<=60);
 subf=f(f<=60);
 subf=fliplr(subf);
 figure();
%  subS=subS./repmat(subf,size(subS,1),1); % whiten
 subS=subS./repmat(max(subS,[],1),size(subS,1),1); % normalize
 [xx yy]=size(subS);
 yy=1:yy;
 xx=1:xx;
 [xi yi]=meshgrid(1:0.1:max(xx),1:0.1:max(yy));
 dataInt=interp2(xx,yy,subS',xi,yi);
 dataInt=dataInt';
 imagesc(flipud(dataInt'));
 for i=1:length(t)
     tstring{i}=t(i);
     tstring{i}=tstring{i}*100;
     tstring{i}=floor(tstring{i});
     tstring{i}=tstring{i}./100;
     tstring{i}=num2str(tstring{i});
 end
 for i=1:length(subf)
     fstring{i}=num2str(floor(subf(i)));
 end
 set(gca,'XTick',1:5*10:size(dataInt',2));
 set(gca,'XTickLabel',tstring(1:5:end));
 set(gca,'YTick',1:3*10:size(dataInt',1));
 set(gca,'YTickLabel',fstring(1:3:end));
 
end