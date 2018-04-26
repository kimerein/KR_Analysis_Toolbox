function makeBurstsByBrainState(outdir,spikes,psth,noThetaTrials)

[fractionpsth,allspikespsth,burstpsth,unitISIplots]=makeBurstPSTH(spikes,psth);
save([outdir '\fractionpsth.mat'],'fractionpsth');
save([outdir '\allspikespsth.mat'],'allspikespsth');
save([outdir '\burstpsth.mat'],'burstpsth');
save([outdir '\unitISIplots.mat'],'unitISIplots');

tri=psth.unitTrials{1};
if length(tri)~=length(noThetaTrials)
    disp('Error');
    return
end
taketri=tri(noThetaTrials);
noTheta_fractionpsth=filtPSTH(fractionpsth,taketri);
save([outdir '\noTheta_fractionpsth.mat'],'noTheta_fractionpsth');
noTheta_allspikespsth=filtPSTH(allspikespsth,taketri);
save([outdir '\noTheta_allspikespsth.mat'],'noTheta_allspikespsth');
noTheta_burstpsth=filtPSTH(burstpsth,taketri);
save([outdir '\noTheta_burstpsth.mat'],'noTheta_burstpsth');

taketri=tri(~noThetaTrials);
theta_fractionpsth=filtPSTH(fractionpsth,taketri);
save([outdir '\theta_fractionpsth.mat'],'theta_fractionpsth');
theta_allspikespsth=filtPSTH(allspikespsth,taketri);
save([outdir '\theta_allspikespsth.mat'],'theta_allspikespsth');
theta_burstpsth=filtPSTH(burstpsth,taketri);
save([outdir '\theta_burstpsth.mat'],'theta_burstpsth');


