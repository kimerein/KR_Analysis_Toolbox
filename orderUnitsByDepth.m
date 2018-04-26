function [order,newassignsinfo2]=orderUnitsByDepth(assignsinfo,trodeOrderByDepth)

% trodeOrderByDepth=[1 2 3 4]; % superficial to deep

[s,orig_i]=sort(assignsinfo.trode);
newassignsinfo.orig_i=orig_i;
newassignsinfo.trode=assignsinfo.trode(orig_i);
newassignsinfo.original_assigns=assignsinfo.original_assigns(orig_i);
newassignsinfo.event_channel=assignsinfo.event_channel(orig_i,:);
for i=1:length(newassignsinfo.trode)
    newassignsinfo.depthTrode(i)=find(trodeOrderByDepth==newassignsinfo.trode(i));
end    

calibrated_evCh=zeros(1,length(assignsinfo.trode));
for i=1:length(assignsinfo.trode)
    curr=assignsinfo.event_channel(i,:);
    [m,ind]=max(curr);
    newcurr=curr;
    newcurr(ind)=-10;
    [m2,ind2]=max(newcurr);
    newcurr2=newcurr;
    newcurr2(ind2)=-10;
    [m3,ind3]=max(newcurr2);
    if curr(ind)==curr(ind2)
        calibrated_evCh(i)=(newassignsinfo.depthTrode(i)-1)*4+mean([ind ind2]);
    elseif curr(ind2)==curr(ind3)
        calibrated_evCh(i)=(newassignsinfo.depthTrode(i)-1)*4+ind;
    else
        a=curr(ind);
        b=curr(ind2);
        y=(100*b)/(a+b);
        x=100-y;
        if ind>ind2
            calibrated_evCh(i)=(newassignsinfo.depthTrode(i)-1)*4+ind-(y/100);
        else
            calibrated_evCh(i)=(newassignsinfo.depthTrode(i)-1)*4+ind+(y/100);
        end
    end
end
[newcal,new_i]=sort(calibrated_evCh);
order=orig_i(new_i);
newassignsinfo2.trode=newassignsinfo.trode(new_i);
newassignsinfo2.original_assigns=newassignsinfo.original_assigns(new_i);
newassignsinfo2.event_channel=newassignsinfo.event_channel(new_i,:);
newassignsinfo2.depthTrode=newassignsinfo.depthTrode(new_i);
newassignsinfo2.calibrated_evCh=newcal;

% figure(); scatter(newassignsinfo.calibrated_evCh,noLedR(order),[],'k'); hold on; scatter(newassignsinfo.calibrated_evCh,ledR(order),[],'r');
% plot([newassignsinfo.calibrated_evCh' newassignsinfo.calibrated_evCh']',[noLedR(order) ledR(order)]','Color',[0.1 0.1 0.1]);