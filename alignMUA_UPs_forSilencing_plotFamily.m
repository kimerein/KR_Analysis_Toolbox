function alignMUA_UPs_forSilencing_plotFamily(expt,mua_xpoints,mua_ypoints,UPstates,useSpontUPs,fileInds,powerRatio)
useLED=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInds));
useStimcond=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInds));

UPthresh=4.6;
stimDuration=5;
ledat=1.3;

runningUPav_noLED=zeros(length(useLED)*20,length(mua_xpoints));
runningUPav_LED=zeros(length(useLED)*20,length(mua_xpoints));
runBefore=zeros(length(useLED)*20,length(mua_xpoints));
runDuring=zeros(length(useLED)*20,length(mua_xpoints));
runDuringnoLED=zeros(1,length(mua_xpoints));
c_before=1;
c_during=1;
c_noLED=1;
c_LED=1;
c_duringNoLED=0;
UPstate_mua=[];
UPstate_pr=[];
ledlocs=[];
pr_xpoints=0:stimDuration/length(powerRatio{1}):stimDuration-(stimDuration/length(powerRatio{1}));
if useSpontUPs==1
    for i=1:length(useLED)
        if useStimcond(i)==9 & ismember(useLED(i),[1 3 5 7])
%         if useStimcond(i)==9 & ismember(useLED(i),[1 3 5 7 4 8])
%         if useStimcond(i)==9 & ismember(useLED(i),[2 4 6 8])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            currpr=powerRatio{i};
            for j=1:size(ups,1)
                if ups(j,1)>1.3 & ups(j,1)<2.1 & ismember(useLED(i),[3 5])
                    continue
                elseif ups(j,1)<1.3 & ups(j,2)>1.25 & ismember(useLED(i),[3 5])
%                 elseif ups(j,1)<1.3 & ups(j,2)>1.25 & ismember(useLED(i),[4 8])
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                        ledlocs=[ledlocs; length(zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))))+find(mua_xpoints<ledat,1,'last')];
                        curr_runBefore=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<ledat)];
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                        ledlocs=[ledlocs; find(mua_xpoints(mua_xpoints>(ups(j,1)-0.5))<ledat,1,'last')];
                        curr_runBefore=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<ledat);
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    end
                    if length(curr_runningUP)>size(runningUPav_LED,2)
                        curr_runningUP=curr_runningUP(1:size(runningUPav_LED,2));
                    elseif length(curr_runningUP)<size(runningUPav_LED,2)
                        curr_runningUP=[curr_runningUP zeros(1,size(runningUPav_LED,2)-length(curr_runningUP))];
                    end
                    if length(curr_runBefore)>size(runningUPav_LED,2)
                        curr_runBefore=curr_runBefore(1:size(runningUPav_LED,2));
                    elseif length(curr_runBefore)<size(runningUPav_LED,2)
                        curr_runBefore=[zeros(1,size(runningUPav_LED,2)-length(curr_runBefore)) curr_runBefore];
                    end
                    if length(curr_runDuring)>size(runningUPav_LED,2)
                        curr_runDuring=curr_runDuring(1:size(runningUPav_LED,2));
                    elseif length(curr_runDuring)<size(runningUPav_LED,2)
                        curr_runDuring=[curr_runDuring zeros(1,size(runningUPav_LED,2)-length(curr_runDuring))];
                    end
                    runningUPav_LED(c_LED,:)=curr_runningUP;
                    runBefore(c_before,:)=curr_runBefore;
                    runDuring(c_during,:)=curr_runDuring;
                    c_LED=c_LED+1;
                    c_before=c_before+1;
                    c_during=c_during+1;
%                 else
                elseif ismember(useLED(i),[1 3 5 7])
%                 elseif ismember(useLED(i),[2 4 6 8])
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                    end
                    if length(curr_runningUP)>size(runningUPav_LED,2)
                        curr_runningUP=curr_runningUP(1:size(runningUPav_LED,2));
                    elseif length(curr_runningUP)<size(runningUPav_LED,2)
                        curr_runningUP=[curr_runningUP zeros(1,size(runningUPav_LED,2)-length(curr_runningUP))];
                    end
                    runningUPav_noLED(c_noLED,:)=curr_runningUP;
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
else
    for i=1:length(useLED)
        if ismember(useStimcond(i),1:8) & ismember(useLED(i),[1 7])
%          if ismember(useStimcond(i),1:8) & ismember(useLED(i),[2 6])
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
                    if length(curr_runningUP)>size(runningUPav_LED,2)
                        curr_runningUP=curr_runningUP(1:size(runningUPav_LED,2));
                    elseif length(curr_runningUP)<size(runningUPav_LED,2)
                        curr_runningUP=[curr_runningUP zeros(1,size(runningUPav_LED,2)-length(curr_runningUP))];
                    end
                    if length(curr_runDuringnoLED)>size(runningUPav_LED,2)
                        curr_runDuringnoLED=curr_runDuringnoLED(1:size(runningUPav_LED,2));
                    elseif length(curr_runDuringnoLED)<size(runningUPav_LED,2)
                        curr_runDuringnoLED=[curr_runDuringnoLED zeros(1,size(runningUPav_LED,2)-length(curr_runDuringnoLED))];
                    end
                    runningUPav_noLED(c_noLED,:)=curr_runningUP;
                    runDuringnoLED=runDuringnoLED+curr_runDuringnoLED;
                    c_noLED=c_noLED+1;
                    c_duringNoLED=c_duringNoLED+1;
                end
            end
%         elseif ismember(useStimcond(i),1:8) & ismember(useLED(i),[3 5])
        elseif ismember(useStimcond(i),1:8) & ismember(useLED(i),[3 5])
%         elseif ismember(useStimcond(i),1:8) & ismember(useLED(i),[4 8])
            ups=UPstates{i};
            currmua_ypoints=mua_ypoints{i};
            for j=1:size(ups,1)
                if ups(j,1)>1 & ups(j,1)<1.3 & ups(j,2)>1.25
                    if isempty(currmua_ypoints(mua_xpoints>=ups(j,1) & mua_xpoints<=ups(j,2)))
                        continue
                    end
                    if ups(j,1)-0.5<0
                        curr_runningUP=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<(ups(j,2)+0.5))];
                        ledlocs=[ledlocs; length(zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))))+find(mua_xpoints<ledat,1,'last')]; 
                        curr_runBefore=[zeros(1,floor(abs(ups(j,1)-0.5)/(mua_xpoints(2)-mua_xpoints(1)))) currmua_ypoints(mua_xpoints>0 & mua_xpoints<ledat)];
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    else
                        curr_runningUP=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<(ups(j,2)+0.5));
                        ledlocs=[ledlocs; find(mua_xpoints(mua_xpoints>(ups(j,1)-0.5))<ledat,1,'last')];
                        curr_runBefore=currmua_ypoints(mua_xpoints>(ups(j,1)-0.5) & mua_xpoints<ledat);
                        curr_runDuring=currmua_ypoints(mua_xpoints>=ledat & mua_xpoints<ledat+3.5);
                    end
                    if length(curr_runningUP)>size(runningUPav_LED,2)
                        curr_runningUP=curr_runningUP(1:size(runningUPav_LED,2));
                    elseif length(curr_runningUP)<size(runningUPav_LED,2)
                        curr_runningUP=[curr_runningUP zeros(1,size(runningUPav_LED,2)-length(curr_runningUP))];
                    end
                    if length(curr_runBefore)>size(runningUPav_LED,2)
                        curr_runBefore=curr_runBefore(1:size(runningUPav_LED,2));
                    elseif length(curr_runBefore)<size(runningUPav_LED,2)
                        curr_runBefore=[zeros(1,size(runningUPav_LED,2)-length(curr_runBefore)) curr_runBefore];
                    end
                    if length(curr_runDuring)>size(runningUPav_LED,2)
                        curr_runDuring=curr_runDuring(1:size(runningUPav_LED,2));
                    elseif length(curr_runDuring)<size(runningUPav_LED,2)
                        curr_runDuring=[curr_runDuring zeros(1,size(runningUPav_LED,2)-length(curr_runDuring))];
                    end
                    runningUPav_LED(c_LED,:)=curr_runningUP;
                    runBefore(c_before,:)=curr_runBefore;
                    runDuring(c_during,:)=curr_runDuring;
                    c_LED=c_LED+1;
                    c_before=c_before+1;
                    c_during=c_during+1;
                end
            end
        end             
    end
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
for i=1:10:c_noLED-1
%     plot(mua_xpoints,smooth(runningUPav_noLED(i,:),10),'Color','k');
    hold on;
end
for i=2:5:c_LED-1
    s=smooth(runningUPav_LED(i,:),10);
    plot(mua_xpoints,s);
    hold on;
    scatter(mua_xpoints(ledlocs(i)),s(ledlocs(i)));
    hold all;
end

figure(); 
plot(mua_xpoints+ledat,smooth(runDuringnoLED,10),'Color','k');
hold on;
for i=1:2:c_before-1
    plot(mua_xpoints,smooth(runBefore(i,:),10));
    hold all;
end
figure(); 
for i=1:2:c_during-1
    plot(mua_xpoints+ledat,smooth(runDuring(i,:),10));
    hold all;
end
% plot(mua_xpoints,smooth(runBefore,10),'Color','g');
% plot(mua_xpoints+ledat,smooth(runDuring,10),'Color','g');