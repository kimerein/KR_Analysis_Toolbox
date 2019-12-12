function concat_all_PSTHs(datadir)

% use averageUnitPSTHs_acrossExpts.m instead

if iscell(datadir)
    d=datadir{1};
    a=load([d '\' 'dLGNpsth']);
    ongoing_psth=a.dLGNpsth;
    for i=2:length(datadir)
        d=datadir{i};
        a=load([d '\' 'dLGNpsth']);
        psth2=a.dLGNpsth;
        ongoing_psth=concatPSTHs(ongoing_psth,psth2);
    end
end


