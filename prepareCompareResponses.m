ledForSweeps=[ledForSweeps1(1:end-1) ledForSweeps(1:end-1)];
LFPbySweep=[LFPbySweep1(1:end-1,:); LFPbySweep(1:end-1,:)];

save('E:\Results\GammaLFPs\KR_2010-09-08\Stim10done\PhotoAligned_Data\all_ledForSweeps.mat','ledForSweeps');
save('E:\Results\GammaLFPs\KR_2010-09-08\Stim10done\PhotoAligned_Data\all_ledAv.mat','ledAv');
save('E:\Results\GammaLFPs\KR_2010-09-08\Stim10done\PhotoAligned_Data\all_LFPbySweep.mat','LFPbySweep');