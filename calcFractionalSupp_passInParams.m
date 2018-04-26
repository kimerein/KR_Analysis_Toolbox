function fracSupp=calcFractionalSupp_passInParams(x,y1,y2,params,spont)

base=params.baseline;
ledWindow=params.ledWindow;
wait=params.wait;
% ledWindow(1)=ledWindow(1)+wait;

if spont==0
    baseActivity1=mean(y1(x>=base(1) & x<=base(2)));
    noLedActivity=mean(y1(x>=ledWindow(1) & x<=ledWindow(2)));
    baseActivity2=mean(y2(x>=base(1) & x<=base(2)));
    ledActivity=mean(y2(x>=ledWindow(1) & x<=ledWindow(2)));
    
    fracSupp=1-((ledActivity-baseActivity2)/(noLedActivity-baseActivity1));
else
    % Don't use baseline
    noLedActivity=mean(y1(x>=ledWindow(1) & x<=ledWindow(2)));
    ledActivity=mean(y2(x>=ledWindow(1) & x<=ledWindow(2)));
    
    fracSupp=1-(ledActivity/noLedActivity);
end