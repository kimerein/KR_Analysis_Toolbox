function [freqs,allp,sum_responses]=getAllFreqResponse_matrix(spikes)

usePeakResponseAmps=0;

a=unique(spikes.assigns);

sum_responses=[];
if usePeakResponseAmps==1
    for i=1:length(a)
        disp(i);
        [freqs,p]=getUnitPeakResponseAmps(spikes,a(i));
        allp(i,:)=p;
    end
else
    for i=1:length(a)
        disp(i);
        [freqs,p,curr_response]=getUnitFreqResponse_matrix(spikes,a(i));
        if i==1
            allp=p;
        else
            allp=allp+p;
        end
        if i==1
            sum_responses=curr_response;
        else
            sum_responses=sum_responses+curr_response;
        end
    end
    sum_responses=sum_responses./length(a);
end

allp(isnan(allp))=0;
figure();
% errorbar(freqs,mean(allp,1),std(allp,[],1));
imagesc(allp);

if usePeakResponseAmps==0
    figure(); 
    imagesc(sum_responses);
end