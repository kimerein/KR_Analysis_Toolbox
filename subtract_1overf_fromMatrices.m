function [newm,newdiags]=subtract_1overf_fromMatrices(m)

doAmp=1;
subtract3xF0=0;
dontSubtractPink=0;
subtractPinkFromAv=0;
peakNorm=0;

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

if subtract3xF0==1
    for i=1:size(m,3)
        m(3,5,i)=m(3,5,i)-m(3,7,i);
        m(3,7,i)=m(3,7,i)-m(3,7,i);
        m(4,7,i)=m(4,7,i)-m(4,10,i);
        m(4,10,i)=m(4,10,i)-m(4,10,i);
        m(11,13,i)=m(11,13,i)-1.5*m(11,15,i);
        m(11,15,i)=m(11,15,i)-m(11,15,i);
        m(12,15,i)=m(12,15,i)-0.5*m(12,12,i);
        m(:,:,i)=m(:,:,i)-min(min(m(:,:,i)));
    end
end

if dontSubtractPink==1
    if doAmp==1
        newm=sqrt(m);
    else
        newm=m;
    end
elseif subtractPinkFromAv==0
    if doAmp==1
        m=sqrt(m);
    end
    newm=nan(size(m));
    for i=1:size(m,3)
%         overf_scale=(m(14,1,i)-mean(m(15,11:15,i)))*(1./sqrt(freqs));
        overf_scale=(mean(m(14:15,1,i))-mean(mean(m(14:15,11:15,i))))*(1./sqrt(freqs));
        newm(:,:,i)=(m(:,:,i)-mean(mean(m(14:15,11:15,i))))-repmat(overf_scale,15,1);
%         figure(); 
%         imagesc(newm(:,:,i));
%         disp('stop');
    end
else
    if doAmp==1
        m=sqrt(m);
    end
    avm=nanmean(m,3);
    overf_scale=(avm(15,1)-mean(avm(15,11:15)))*(1./sqrt(freqs));
    newm=nan(size(m));
    for i=1:size(m,3)
        newm(:,:,i)=m(:,:,i);
        newm(:,:,i)=(newm(:,:,i)-mean(avm(15,11:15)))-repmat(overf_scale,15,1);
    end
end

newdiags=nan(size(m,3),15);
for i=1:size(m,3)
    for j=1:15
        newdiags(i,j)=newm(j,j,i);
    end
end

if peakNorm==1
    for i=1:size(newdiags,1)
        newdiags(i,:)=newdiags(i,:)./max(newdiags(i,:));
        newm(:,:,i)=newm(:,:,i)./max(newdiags(i,:));
    end
end

