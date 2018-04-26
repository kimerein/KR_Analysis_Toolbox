function [A_vector,tau_vector,normOfResid_vector,bestWaitS]=acrossLayersTau_wrapper(expt,useTheseTrodes,trodeChs,layerSpecs)

waitS=0:0.001:0.01;
for i=1:length(waitS)
    disp(i);
    [A{i},tau{i},normOfResid{i}]=calculate_tau_across_layers(expt,[],[],[],waitS(i));
end

for i=1:length(waitS)
    a=tau{i};
%     b=[a{1} a{2} a{3} a{4} a{5}];
    b=a{1};
    tau{i}=b;
    a=A{i};
%     b=[a{1} a{2} a{3} a{4} a{5}];
    b=a{1};
    A{i}=b;
    a=normOfResid{i};
%     b=[a{1} a{2} a{3} a{4} a{5}];
    b=a{1};
    normOfResid{i}=b;
end

% Don't want algorithm to choose taus that are infinity
% so set A/NOR as 0 if tau is infinity
for i=1:length(waitS)
%     tauinf=isinf(tau{i});
    tauinf=tau{i}<0;
    temp=A{i};
    temp(tauinf)=0;
    A{i}=temp;
end

greatestAoverNOR=max(A{1}./normOfResid{1},A{2}./normOfResid{2});
firstGreaterThanNext=A{1}./normOfResid{1}>=A{2}./normOfResid{2};
bestWaitS=zeros(size(A{1}));
bestWaitS(firstGreaterThanNext)=1;
bestWaitS(~firstGreaterThanNext)=2;
for i=3:length(waitS)
    nextGreaterThanFirst=greatestAoverNOR<(A{i}./normOfResid{i});
    greatestAoverNOR=max(greatestAoverNOR,A{i}./normOfResid{i});
    bestWaitS(nextGreaterThanFirst)=i;
end

A_vector=zeros(size(A{1}));
tau_vector=zeros(size(A{1}));
normOfResid_vector=zeros(size(A{1}));
for i=1:length(waitS)
    currA=A{i};
    currTau=tau{i};
    currNOR=normOfResid{i};
    A_vector(bestWaitS==i)=currA(bestWaitS==i);
    tau_vector(bestWaitS==i)=currTau(bestWaitS==i);
    normOfResid_vector(bestWaitS==i)=currNOR(bestWaitS==i);
end


    
    