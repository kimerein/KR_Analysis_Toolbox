function powerInBand=takeFrequencyBand(specgram,freqBand,noTheta)

if length(size(specgram))>2
    powerInBand=nan(size(specgram,3),size(specgram,1)); 
    for i=1:size(specgram,3)
        powerInBand(i,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2),i),2),size(specgram,1),1);
    end
else
    powerInBand=nan(1,size(specgram,1)); 
    powerInBand(1,:)=reshape(nanmean(specgram(:,noTheta.allS.f>=freqBand(1) & noTheta.allS.f<=freqBand(2)),2),size(specgram,1),1);
end

end