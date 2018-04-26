function [thetaDiff,S,t,f]=plotHippoTheta(LFPbySweep,ledConds,stimConds,usel,uses,Fs,S,t,f,thetaDiff)

params.tapers=[3 5];
params.Fs=Fs;
params.fpass=[0.5 2000];
params.trialave=0;
movingwin=[0.8 0.2];

thetaLowband=[2 4];
thetaHighband=[5 7];
norm=1;

if isempty(S)
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
else
end

[n,xout]=hist(thetaDiff(1:end),300);
figure(); 
plot(xout,n);

% Plot example trial
% trialn=10:size(S,3);
% trialn=85:105;
trialn=91:94;
disp(size(S,3));
concatS=[];
concatThetaDiff=[];
concatT=[];
k=1;
for j=1:length(trialn)
    i=trialn(j);
    concatS=[concatS reshape(S(:,:,i),size(S,1),size(S,2))'];
    concatThetaDiff=[concatThetaDiff thetaDiff(i,:)];
    concatT=[concatT t+(k-1)*max(t)];
    k=k+1;
end
ds=1;
figure();
imagesc(downSampAv(concatT,ds),f,downSampMatrix(concatS,ds));
ylim([min(f) 30]);
colorbar();
figure();
plot(downSampAv(concatT,ds),downSampAv(concatThetaDiff,ds));
colorbar();

% figure();
% outim.t=t;
% outim.f=f;
% outim.im=reshape(S(:,:,trialn),size(S,1),size(S,2))';
% imagesc(t,f,outim.im);
% 
% figure();
% plot(t,thetaDiff(trialn,:));

