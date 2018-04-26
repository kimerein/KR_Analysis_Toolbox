function output=cycleAv_wrapper(expt,spikes,nBins,output,showFreqInd)

peakAlign=1;

a=unique(spikes.assigns);
% useTheseAssigns=find(ismember(a,[2 4 6 9 10 11 12 25 26 27 28 29 30 31 37 38 39 40 41 42 43 45 48 51 56 57 59 62 63 85 86 88 102 106 108 114 115 116 117 118 119 120 127]));
useTheseAssigns=1:length(a);
stimFreqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
amberFreqs=[1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];

if isempty(output)
%     for i=1:length(a)
    for i=useTheseAssigns
        disp(i);
        [cav_con,cav_led,sOverBas_con,sOverBas_led]=plotCycleAverage(expt,spikes,a(i),nBins,1);
        all_sOverBas_con(i,:)=sOverBas_con;
        all_sOverBas_led(i,:)=sOverBas_led;
        for j=1:length(cav_con)
            ccav_con=cav_con{j};
            all_con_x(j).unitCycs(i,:)=ccav_con(1,:);
            all_con_y(j).unitCycs(i,:)=ccav_con(2,:);
            ccav_led=cav_led{j};
            all_led_x(j).unitCycs(i,:)=ccav_led(1,:);
            all_led_y(j).unitCycs(i,:)=ccav_led(2,:);
        end
    end  
    for i=1:length(stimFreqs)
        all_con_x(i).m=mean(all_con_x(i).unitCycs,1);
        all_con_x(i).s=std(all_con_x(i).unitCycs,[],1);
        all_con_y(i).m=mean(all_con_y(i).unitCycs,1);
        all_con_y(i).s=std(all_con_y(i).unitCycs,[],1);
        all_led_x(i).m=mean(all_led_x(i).unitCycs,1);
        all_led_x(i).s=std(all_led_x(i).unitCycs,[],1);
        all_led_y(i).m=mean(all_led_y(i).unitCycs,1);
        all_led_y(i).s=std(all_led_y(i).unitCycs,[],1);
    end
    output.all_con_x=all_con_x;
    output.all_con_y=all_con_y;
    output.all_led_x=all_led_x;
    output.all_led_y=all_led_y;
else
    all_con_x=output.all_con_x;
    all_con_y=output.all_con_y;
    all_led_x=output.all_led_x;
    all_led_y=output.all_led_y;
end

if ~isempty(useTheseAssigns)
    for i=1:length(all_con_x)
        all_con_x(i).unitCycs=all_con_x(i).unitCycs(useTheseAssigns,:);
        all_con_y(i).unitCycs=all_con_y(i).unitCycs(useTheseAssigns,:);
        all_led_x(i).unitCycs=all_led_x(i).unitCycs(useTheseAssigns,:);
        all_led_y(i).unitCycs=all_led_y(i).unitCycs(useTheseAssigns,:);
    end
    figure(); 
    errorbar(stimFreqs,mean(all_sOverBas_con,1),std(all_sOverBas_con,[],1),'Color','k');
    hold on;
    errorbar(stimFreqs,mean(all_sOverBas_led,1),std(all_sOverBas_led,[],1),'Color','r');
    title('Average Unit Freq. Resp. -- mean during stim. over mean during base.');
end

mid=5;
if peakAlign==1
    % Peak align cycle averages
    for i=1:length(stimFreqs)
        n=floor(size(all_con_y(i).unitCycs,2)/2);
        for j=1:size(all_con_y(i).unitCycs,1)
            base=all_con_y(i).unitCycs(j,1:n);
            stim=all_con_y(i).unitCycs(j,n+1:end);
            [basePeak,basePeakInd]=max(base);
            [stimPeak,stimPeakInd]=max(stim);
            newbase=[base(basePeakInd:end) base(1:basePeakInd-1)];
            newbase=[newbase(mid+1:end) newbase(1:mid)];
            newstim=[stim(stimPeakInd:end) stim(1:stimPeakInd-1)];
            newstim=[newstim(mid+1:end) newstim(1:mid)];
            all_con_y(i).unitCycs(j,1:n)=newbase;
            all_con_y(i).unitCycs(j,n+1:end)=newstim;
            base=all_led_y(i).unitCycs(j,1:n);
            stim=all_led_y(i).unitCycs(j,n+1:end);
            newbase=[base(basePeakInd:end) base(1:basePeakInd-1)];
            newbase=[newbase(mid+1:end) newbase(1:mid)];
            newstim=[stim(stimPeakInd:end) stim(1:stimPeakInd-1)];
            newstim=[newstim(mid+1:end) newstim(1:mid)];
            all_led_y(i).unitCycs(j,1:n)=newbase;
            all_led_y(i).unitCycs(j,n+1:end)=newstim;
        end
    end
    for i=1:length(stimFreqs)
        all_con_x(i).m=mean(all_con_x(i).unitCycs,1);
        all_con_x(i).s=std(all_con_x(i).unitCycs,[],1);
        all_con_y(i).m=mean(all_con_y(i).unitCycs,1);
        all_con_y(i).s=std(all_con_y(i).unitCycs,[],1);
        all_led_x(i).m=mean(all_led_x(i).unitCycs,1);
        all_led_x(i).s=std(all_led_x(i).unitCycs,[],1);
        all_led_y(i).m=mean(all_led_y(i).unitCycs,1);
        all_led_y(i).s=std(all_led_y(i).unitCycs,[],1);
    end
    output.all_con_x=all_con_x;
    output.all_con_y=all_con_y;
    output.all_led_x=all_led_x;
    output.all_led_y=all_led_y;
end
            
figure(); 
i=showFreqInd;
errorbar(all_con_x(i).m,all_con_y(i).m,all_con_y(i).s,'Color','k');
hold on; 
errorbar(all_led_x(i).m,all_led_y(i).m,all_led_y(i).s,'Color','r');