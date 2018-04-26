l1=[ledAmp1(1,:); ledAmp2(1,:); ledAmp3(1,:); ledAmp4(1,:)];
l2=[ledAmp1(2,:); ledAmp2(2,:); ledAmp3(2,:); ledAmp4(2,:)];
l3=[ledAmp1(3,:); ledAmp2(3,:); ledAmp3(3,:); ledAmp4(3,:)];
l4=[ledAmp1(4,:); ledAmp2(4,:); ledAmp3(4,:); ledAmp4(4,:)];

l1=l1-[l1(1,:)+1; l1(1,:); l1(1,:); l1(1,:)];
l2=l2-[l2(1,:)+1; l2(1,:); l2(1,:); l2(1,:)];
l3=l3-[l3(1,:)+1; l3(1,:); l3(1,:); l3(1,:)];
l4=l4-[l4(1,:)+1; l4(1,:); l4(1,:); l4(1,:)];

figure();
plot(l1);
hold on;
plot(l2);
plot(l3);
plot(l4);