function showSpontaneousDaq(daqFileName)

% Add option to also dispay photodiode channel!!!

persistent data time photoData ledData currDaqFile

physChannel=11;
photoChannel=6;
ledChannel=4;
chunkLengthInSecs=50; 
chunkN=1;

showPhotoLED=1;

if isempty(data) || isempty(time) || isempty(currDaqFile) || ~strcmp(currDaqFile,daqFileName)
    [data, time]=daqread(daqFileName,'Triggers',[1 16]);
    %[data, time]=daqread(daqFileName);
    thisData=data(:,physChannel);
    photoData=data(:,photoChannel);
    ledData=data(:,ledChannel);
    data=thisData;
    clear thisData;
    currDaqFile=daqFileName;
end

chunkLengthInInds=find(floor(time)==chunkLengthInSecs,1,'first');

data=data(~isnan(data));
photoData=photoData(~isnan(photoData));
ledData=ledData(~isnan(ledData));
time=time(~isnan(time));

for i=chunkN
    scrsz=get(0,'ScreenSize');
    figure('Position',[100 scrsz(4)/2-300 scrsz(3)*0.75 scrsz(4)/6]);
    %plot(time(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i),data(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i),'Color',[rand(1) rand(1) rand(1)]);
    plot(time(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i),data(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i)-max(data(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i)),'Color','black'); % Note times 0.5 is KR hack
    if showPhotoLED==1
        hold on;
        photoRange=max(photoData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i))-min(photoData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i));
        sigRange=abs(max(data(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i))-min(data(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i)));
        a=(1/10)*(sigRange/photoRange);
        photoSignal2=a*photoData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i)-min(a*photoData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i));
        plot(time(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i),photoSignal2,'Color',[0.1 0.1 0.1]);
        ledRange=max(ledData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i))-min(ledData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i));
        a=(1/10)*(sigRange/ledRange);
        ledSignal2=a*ledData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i)-min(a*ledData(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i));
        plot(time(chunkLengthInInds*(i-1)+1:chunkLengthInInds*i),ledSignal2,'Color','cyan');
    end
    axis([time(chunkLengthInInds*(i-1)+1) time(chunkLengthInInds*i) 0 1]);
    axis 'auto y';
end
