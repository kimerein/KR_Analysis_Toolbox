
spontSeries(1,1:4)=[fSpont(2:4); fSpont(1)];
spontSeries(2,1:4)=[fSpont(6:8); fSpont(5)];
spontSeries(3,1:4)=[fSpont(10:12); fSpont(9)];
spontSeries(4,1:4)=[fSpont(13:15); nan];
spontSeries(5,1:4)=[fSpont(17:19); fSpont(16)];
spontSeries(6,1:4)=[fSpont(21:23); fSpont(20)];
spontSeries(7,1:4)=[fSpont(25:27); fSpont(24)];


evSeries(1,1:4)=[fEv(2:4); fEv(1)];
evSeries(2,1:4)=[fEv(6:8); fEv(5)];
evSeries(3,1:4)=[fEv(10:12); fEv(9)];
evSeries(4,1:4)=[fEv(13:15); nan];
evSeries(5,1:4)=[fEv(17:19) ;fEv(16)];
evSeries(6,1:4)=[fEv(21:23); fEv(20)];
evSeries(7,1:4)=[fEv(25:27) ;fEv(24)];

% evSeries(1,1:4)=fEv(1:4);
% evSeries(2,1:4)=fEv(5:8);
% evSeries(3,1:4)=fEv(9:12);
% evSeries(4,1:4)=[fEv(13:15); nan];
% evSeries(5,1:4)=fEv(16:19);
% evSeries(6,1:4)=fEv(20:23);
% evSeries(7,1:4)=fEv(24:27);

figure();
scatter(1-evSeries(1:end-2,:),1-spontSeries(1:end-2,:),[0 0 0]);
hold on;

% for i=1:size(spontSeries-2,1)
%     plot(1-evSeries(i,:),1-spontSeries(i,:),'Color','k');
%     hold on;
% end
% for i=1:2
%     plot(1-evSeries(5+i,:),1-spontSeries(5+i,:),'Color','red');
%     hold on;
% end