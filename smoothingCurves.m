% for i=1:length(Ydata1)-1
%     newYdata1(i)=mean(Ydata1(i:i+1));
% end
    
a1=n2Ydata1;
a2=n2Ydata2;
x1=n2Xdata;

for i=1:length(a1)-1
    b1(i)=mean(a1(i:i+1));
end
for i=1:length(a2)-1
    b2(i)=mean(a2(i:i+1));
end
for i=1:length(x1)-1
    c(i)=mean(x1(i:i+1));
end
figure(); plot(c,b1,'Color','k'); hold on; plot(c,b2,'Color','r');

n2primeYdata1=b1;
n2primeYdata2=b2;
n2primeXdata=c;