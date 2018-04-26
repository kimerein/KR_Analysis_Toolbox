function [poppvalmatrix,finalpval]=comparePopVecPvalMatrix(spikes,assigns,ledConds,window)

stimconds=[1:9];

poppvalmatrix=zeros(length(stimconds),length(stimconds));
withinCondDistances=[];
acrossCondDistances=[];

withinLibrary=cell(length(stimconds),1);
for i=1:length(stimconds)
    for j=1:length(stimconds)
        disp(i);
        disp(j);
        if ~isempty(withinLibrary{i})
            useL=withinLibrary{i};
        else
            useL=[];
        end
        if ~isempty(withinLibrary{j})
            useA=withinLibrary{j};
        else
            useA=[];
        end
        [withinCond,acrossCond,pvalMatrix,pvalForPopvecs,pvL,pvC]=comparePopulationVectors(spikes,assigns,i,ledConds,j,ledConds,window,[],2,useL,useA);
        withinLibrary{i}=pvL;
        withinLibrary{j}=pvC;
        poppvalmatrix(i,j)=pvalForPopvecs;
        withinCondDistances=[withinCondDistances; withinCond.distances];
        acrossCondDistances=[acrossCondDistances; acrossCond.distances];
    end
end

[n1,x1]=hist(withinCondDistances,10);
[n2,x2]=hist(acrossCondDistances,10);
figure(); plot(x1,n1,'Color','black'); hold on; plot(x2,n2,'Color','red');
[n1,x1]=histnorm(withinCondDistances,10);
[n2,x2]=histnorm(acrossCondDistances,10);
figure(); plot(x1,n1,'Color','black'); hold on; plot(x2,n2,'Color','red');
disp('final pval');
finalpval=ranksum(withinCondDistances,acrossCondDistances);
disp(finalpval);

