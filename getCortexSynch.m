function cxSynch=getCortexSynch(LFPbySweep,ledConds,stimConds,usel,uses,Fs)

takeOnlyNoLED=1;

params.tapers=[3 5];
params.Fs=Fs;
params.fpass=[0.5 2000];
params.trialave=0;
movingwin=[0.8 0.2];

cxLowband=[0.5 6];
cxHighband=[30 80];
norm=1;

L=LFPbySweep{1};
L=L(ismember(ledConds,usel) & ismember(stimConds,uses),:);
[S,t,f]=mtspecgramc(L',movingwin,params);

if norm==1
    for i=1:size(S,3)
        for j=1:size(S,1)
%             S(j,:,i)=S(j,:,i)./sum(S(j,:,i));
            S(j,:,i)=S(j,:,i)./sum(S(j,:,i));
        end
    end
end

if takeOnlyNoLED==1
    cxSynch=nan(size(S,3),size(S(t<1.2 | t>2.5,:,:),1));
else
    cxSynch=nan(size(S,3),size(S,1));
end
for i=1:size(S,3)
%     if mod(i,10)==0
%         disp(i);
%     end
if takeOnlyNoLED==0
    cxSynch(i,:)=reshape(nanmean(S(:,f>=cxLowband(1) & f<=cxLowband(2),i),2)-nanmean(S(:,f>=cxHighband(1) & f<=cxHighband(2),i),2),1,size(S,1));
else
    cxSynch(i,:)=reshape(nanmean(S(t<1.2 | t>2.5,f>=cxLowband(1) & f<=cxLowband(2),i),2)-nanmean(S(t<1.2 | t>2.5,f>=cxHighband(1) & f<=cxHighband(2),i),2),1,size(S(t<1.2 | t>2.5,:,:),1));
end
end

[n,xout]=hist(cxSynch(1:end),300);
figure(); 
plot(xout,n);

[n,xout]=hist(nanmean(cxSynch,2),100);
figure(); 
plot(xout,n);