function alignedResp=alignFreqResponses(freqs,resp,template)

alignToTemplate=0;
alignToFreqs=[1 2];
% Align to template
if isempty(template)
    template=resp(1,:);
end

% Assume zero is zero
alignedResp=zeros(size(resp));
scalesToTry=0.5:0.01:3;
for i=1:size(resp,1)
    if alignToTemplate==1
        currResp=resp(i,:);
        bestSummedErr=10000;
        bestScale=1;
        for j=1:length(scalesToTry)
            currScale=scalesToTry(j);
            newVers=currResp.*currScale;
            if sum(abs(newVers-template))<bestSummedErr
                bestScale=currScale;
                bestSummedErr=sum(abs(newVers-template));
            end
        end
        alignedResp(i,:)=currResp.*bestScale;
    else
        currResp=resp(i,:);
%         bestScale=1/mean(currResp(ismember(freqs,alignToFreqs)));
%         bestScale=1/sum(currResp);
        bestScale=1/max(currResp);
        alignedResp(i,:)=currResp.*bestScale;
    end
end

figure(); 
errorbar(freqs,mean(alignedResp,1),std(alignedResp,[],1)./sqrt(size(alignedResp,1))); % std error
    