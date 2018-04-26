function alignMUA_UPs_forThalSilencing(expt,mua_xpoints,mua_ypoints,UPstates,useSpontUPs,fileInds,powerRatio)
useLED=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInds));
useStimcond=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInds));

UPthresh=4.6;
stimDuration=5;
ledat=1.3;

runningUPav_noLED=zeros(1,length(mua_xpoints));
runningUPav_LED=zeros(1,length(mua_xpoints));
runBefore=zeros(1,length(mua_xpoints));
runDuring=zeros(1,length(mua_xpoints));
runDuringnoLED=zeros(1,length(mua_xpoints));
c_before=0;
c_during=0;
c_noLED=0;
c_LED=0;
c_duringNoLED=0;
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
                elseif ups(j,1)<1.3 & ups(j,2)>1.25 & ismember(useLED(i),[3 5])
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                        curr_runBefore=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<ledat)];
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                        curr_runBefore=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<ledat);
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    end
                    if length(curr_runningUP)>length(runningUPav_LED)
                        curr_runningUP=curr_runningUP(1:length(runningUPav_LED));
                    elseif length(curr_runningUP)<length(runningUPav_LED)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav_LED)-length(curr_runningUP))];
                    end
                    if length(curr_runBefore)>length(runningUPav_LED)
                        curr_runBefore=curr_runBefore(1:length(runningUPav_LED));
                    elseif length(curr_runBefore)<length(runningUPav_LED)
                        curr_runBefore=[zeros(1,length(runningUPav_LED)-length(curr_runBefore)) curr_runBefore];
                    end
                    if length(curr_runDuring)>length(runningUPav_LED)
                        curr_runDuring=curr_runDuring(1:length(runningUPav_LED));
                    elseif length(curr_runDuring)<length(runningUPav_LED)
                        curr_runDuring=[curr_runDuring zeros(1,length(runningUPav_LED)-length(curr_runDuring))];
                    end
                    runningUPav_LED=runningUPav_LED+curr_runningUP;
                    runBefore=runBefore+curr_runBefore;
                    runDuring=runDuring+curr_runDuring;
                    c_LED=c_LED+1;
                    c_before=c_before+1;
                    c_during=c_during+1;
                else
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                    end
                    if length(curr_runningUP)>length(runningUPav_LED)
                        curr_runningUP=curr_runningUP(1:length(runningUPav_LED));
                    elseif length(curr_runningUP)<length(runningUPav_LED)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav_LED)-length(curr_runningUP))];
                    end
                    runningUPav_noLED=runningUPav_noLED+curr_runningUP;
                    c_noLED=c_noLED+1;
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
    runningUPav_LED=runningUPav_LED./c_LED;
    runningUPav_noLED=runningUPav_noLED./c_noLED;
    runBefore=runBefore./c_before;
    runDuring=runDuring./c_during;
else
    for i=1:length(useLED)
        if ismember(useStimcond(i),1:8) & ismember(useLED(i),[1 7])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            for j=1:size(ups,1)
%                 if ups(j,1)>1 & ups(j,1)<1.5
                if ups(j,1)>1 & ups(j,1)<1.3
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                        curr_runDuringnoLED=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                        curr_runDuringnoLED=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    end
                    if length(curr_runningUP)>length(runningUPav_LED)
                        curr_runningUP=curr_runningUP(1:length(runningUPav_LED));
                    elseif length(curr_runningUP)<length(runningUPav_LED)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav_LED)-length(curr_runningUP))];
                    end
                    if length(curr_runDuringnoLED)>length(runningUPav_LED)
                        curr_runDuringnoLED=curr_runDuringnoLED(1:length(runningUPav_LED));
                    elseif length(curr_runDuringnoLED)<length(runningUPav_LED)
                        curr_runDuringnoLED=[curr_runDuringnoLED zeros(1,length(runningUPav_LED)-length(curr_runDuringnoLED))];
                    end
                    runningUPav_noLED=runningUPav_noLED+curr_runningUP;
                    runDuringnoLED=runDuringnoLED+curr_runDuringnoLED;
                    c_noLED=c_noLED+1;
                    c_duringNoLED=c_duringNoLED+1;
                end
            end
        elseif ismember(useStimcond(i),1:8) & ismember(useLED(i),[3 5])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            for j=1:size(ups,1)
                if ups(j,1)>1 & ups(j,1)<1.3 & ups(j,2)>1.25
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                        curr_runBefore=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<ledat)];
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                        curr_runBefore=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<ledat);
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    end
                    if length(curr_runningUP)>length(runningUPav_LED)
                        curr_runningUP=curr_runningUP(1:length(runningUPav_LED));
                    elseif length(curr_runningUP)<length(runningUPav_LED)
                        curr_runningUP=[curr_runningUP zeros(1,length(runningUPav_LED)-length(curr_runningUP))];
                    end
                    if length(curr_runBefore)>length(runningUPav_LED)
                        curr_runBefore=curr_runBefore(1:length(runningUPav_LED));
                    elseif length(curr_runBefore)<length(runningUPav_LED)
                        curr_runBefore=[zeros(1,length(runningUPav_LED)-length(curr_runBefore)) curr_runBefore];
                    end
                    if length(curr_runDuring)>length(runningUPav_LED)
                        curr_runDuring=curr_runDuring(1:length(runningUPav_LED));
                    elseif length(curr_runDuring)<length(runningUPav_LED)
                        curr_runDuring=[curr_runDuring zeros(1,length(runningUPav_LED)-length(curr_runDuring))];
                    end
                    runningUPav_LED=runningUPav_LED+curr_runningUP;
                    runBefore=runBefore+curr_runBefore;
                    runDuring=runDuring+curr_runDuring;
                    c_LED=c_LED+1;
                    c_before=c_before+1;
                    c_during=c_during+1;
                end
            end
        end             
    end
    runningUPav_LED=runningUPav_LED./c_LED;
    runningUPav_noLED=runningUPav_noLED./c_noLED;
    runBefore=runBefore./c_before;
    runDuring=runDuring./c_during;
    runDuringnoLED=runDuringnoLED./c_duringNoLED;
end

% if useSpontUPs==1
%     figure(); 
%     [n,x]=hist(UPstate_mua,10);
%     plot(x,n);
%     figure(); 
%     [n,x]=hist(UPstate_pr,10);
%     plot(x,n);
%     figure(); 
%     scatter(UPstate_mua,UPstate_pr);
% end

figure(); 
plot(mua_xpoints,smooth(runningUPav_noLED,10),'Color','k');
hold on;
plot(mua_xpoints,smooth(runningUPav_LED,10),'Color','r');
figure(); 
plot(mua_xpoints,smooth(runningUPav_noLED,10),'Color','k');
hold on;
plot(mua_xpoints+ledat,smooth(runDuringnoLED,10),'Color','k');
plot(mua_xpoints,smooth(runBefore,10),'Color','r');
plot(mua_xpoints+ledat,smooth(runDuring,10),'Color','r');