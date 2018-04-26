figure(); 
plot(exptData.xpoints,exptData.ypoints1,'Color','k'); 
hold on; 
plot(exptData.xpoints,exptData.ypoints2,'Color','r'); 

%%
base=[0 1]; 
su=[1.3 1.425]; 
fras=1-(mean(exptData.ypoints2(exptData.xpoints>su(1) & exptData.xpoints<su(2)))-mean(exptData.ypoints2(exptData.xpoints>base(1) & exptData.xpoints<base(2))))/(mean(exptData.ypoints1(exptData.xpoints>su(1) & exptData.xpoints<su(2)))-mean(exptData.ypoints1(exptData.xpoints>base(1) & exptData.xpoints<base(2))));
