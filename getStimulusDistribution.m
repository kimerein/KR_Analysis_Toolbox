function getStimulusDistribution()

offset=100;
x=0:0.0001:10; % in s
s=100*cos(x*2*pi*2)+offset; % for temporal frequency = 2 Hz
% figure(); 
% plot(x,s);

phases=0:0.01:2*pi;
allavs=[];
for k=1:length(phases)
    s=100*cos(x*2*pi*2+phases(k))+offset;
    bins=1:0.25:9;
    savs=zeros(1,length(bins)-1);
    for j=1:length(bins)-1
        savs(j)=mean(s(x>bins(j) & x<bins(j+1)));
    end
    allavs=[allavs savs];
end
    
[n,x]=hist(allavs);
figure(); 
plot(x,n);