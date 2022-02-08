function plotFilledCDF(data)

figure(); 
[f,x,flo,fup]=ecdf(data,'Bounds','on');
temp=find(x>=0,1,'first');
fill([fliplr(x(x>0)); zeros(4,1); x(x>0)],[ones(size(f(x>0))); [1; 0.98; 0.97; f(temp)]; f(x>0)],[0.5 0.5 0.5]);
hold on;
temp=find(x<=0,1,'last');
fill([fliplr(x(x<0)); x(x<0); zeros(4,1)],[zeros(size(f(x<0))); f(x<0);  [f(temp); 0.02; 0.01; 0]],[0.9 0.5 0.5]);
plot(x,f,'Color','k');


end