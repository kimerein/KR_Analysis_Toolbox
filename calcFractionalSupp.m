function fracSupp=calcFractionalSupp(x,y1,y2)

base=[0.6 1];
ledWindow=[3.3 3.8];
wait=0.04;
ledWindow(1)=ledWindow(1)+wait;

baseActivity1=mean(y1(x>=base(1) & x<=base(2)));
noLedActivity=mean(y1(x>=ledWindow(1) & x<=ledWindow(2)));
baseActivity2=mean(y2(x>=base(1) & x<=base(2)));
ledActivity=mean(y2(x>=ledWindow(1) & x<=ledWindow(2)));

fracSupp=1-((ledActivity-baseActivity2)/(noLedActivity-baseActivity1));