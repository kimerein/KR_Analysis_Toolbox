function alignMUA_UPstates(expt,mua_xpoints,mua_ypoints,UPstates,useSpontUPs,fileInds,powerRatio)
useLED=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInds));
useStimcond=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInds));

UPthresh=4.6;
stimDuration=5;

runningUPav=zeros(1,length(mua_xpoints));
c=0;
UPstate_mua=[];
UPstate_pr=[];
pr_xpoints=0:stimDuration/length(powerRatio{1}):stimDuration-(stimDuration/length(powerRatio{1}));
if useSpontUPs==1
    for i=1:length(useLED)
        if useStimcond(i)==9 & ismember(useLED(i),[1 3 5 7])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            currpr=powerRatio{i};
            for j=1:size(ups,1)
                if ups(j,1)>1.3 & ups(j,1)<2.1 & ismember(useLED(i),[3 5])
                    continue
                elseif ups(j,1)<1.3 & ups(j,2)>1.3 & ismember(useLED(i),[3 5])
                    continue
                else
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                    end
                    if length(curr_runningUP)>length(runningUPav)
                        curr_runningUP=curr_runningUP(1:length(runningUPav));
                    elseif length(curr_runningUP)<length(runningUPav)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav)-length(curr_runningUP))];
                    end
                    runningUPav=runningUPav+curr_runningUP;
                    c=c+1;
                    upmua=mean(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)));
                    uppr=mean(currpr(pr_xpoints>=ups(j,1) & pr_xpoints<=ups(j,2)));
                    if isnan(upmua)
                        disp('stop');
                    end
                    if isnan(uppr)
                        disp('stop');
                    end
                    UPstate_mua=[UPstate_mua upmua];
                    UPstate_pr=[UPstate_pr uppr];
                end
            end
        elseif ismember(useStimcond(i),1:8) & ismember(useLED(i),[1 3 5 7])
%             ups=UPstates{i};
%             currmua_ypoints=mua_ypoints{i};
%             for j=1:size(ups,1)
%                 if (ups(j,1)<1 & ups(j,2)<1) || ups(j,1)>4.5
%                     if ups(j,1)-0.5<0
%                         curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
%                     else
%                         curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
%                     end
%                     if length(curr_runningUP)>length(runningUPav)
%                         curr_runningUP=curr_runningUP(1:length(runningUPav));
%                     elseif length(curr_runningUP)<length(runningUPav)
%                         curr_runningUP=[curr_runningUP zeros(1,length(runningUPav)-length(curr_runningUP))];
%                     end
%                     runningUPav=runningUPav+curr_runningUP;
%                     c=c+1;
%                 else
%                     continue
%                 end
%             end
        end
    end
    runningUPav=runningUPav./c;
else
    for i=1:length(useLED)
        if ismember(useStimcond(i),1:8) & ismember(useLED(i),[1 7])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            for j=1:size(ups,1)
                if ups(j,1)>1 & ups(j,1)<1.3
%                 if ups(j,1)<1 
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                    end
                    if length(curr_runningUP)>length(runningUPav)
                        curr_runningUP=curr_runningUP(1:length(runningUPav));
                    elseif length(curr_runningUP)<length(runningUPav)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav)-length(curr_runningUP))];
                    end
                    runningUPav=runningUPav+curr_runningUP;
                    c=c+1;
                end
            end
        end 
    end
    runningUPav=runningUPav./c;
end

if useSpontUPs==1
    figure(); 
    [n,x]=hist(UPstate_mua,10);
    plot(x,n);
    figure(); 
    [n,x]=hist(UPstate_pr,10);
    plot(x,n);
    figure(); 
    scatter(UPstate_mua,UPstate_pr);
end

figure(); 
plot(mua_xpoints,smooth(runningUPav,10));