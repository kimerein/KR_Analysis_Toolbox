function thetaDiff=getHippoTheta(LFPbySweep,ledConds,stimConds,usel,uses,Fs)

params.tapers=[3 5];
params.Fs=Fs;
params.fpass=[0.5 2000];
params.trialave=0;
movingwin=[0.8 0.2];

thetaLowband=[2 4];
thetaHighband=[5 7];
norm=1;

L=LFPbySweep{1};

if any(usel-floor(usel)>0.0001)
else
    L=L(ismember(ledConds,usel) & ismember(stimConds,uses),:);
end
[S,t,f]=mtspecgramc(L',movingwin,params);

if norm==1
    for i=1:size(S,3)
        for j=1:size(S,1)
%             S(j,:,i)=S(j,:,i)./sum(S(j,:,i));
            S(j,:,i)=S(j,:,i)./sum(S(j,:,i));
        end
    end
end

thetaDiff=nan(size(S,3),size(S,1));
for i=1:size(S,3)
%     if mod(i,10)==0
%         disp(i);
%     end
    thetaDiff(i,:)=reshape(nanmean(S(:,f>=thetaHighband(1) & f<=thetaHighband(2),i),2)-nanmean(S(:,f>=thetaLowband(1) & f<=thetaLowband(2),i),2),1,size(S,1));
%     thetaDiff(i,:)=reshape(nanmean(S(:,f>=thetaHighband(1) & f<=thetaHighband(2),i),2),1,size(S,1));
end

[n,xout]=hist(thetaDiff(1:end),300);
figure(); 
plot(xout,n);