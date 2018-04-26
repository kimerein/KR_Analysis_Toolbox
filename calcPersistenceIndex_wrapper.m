function [return_pers_index,a,b,c,d,gain_persIndex,ev,b_over_a,base_over_a]=calcPersistenceIndex_wrapper()

% dataDir{1}='W:\Analysis Computer\Persistent Units\Mawake54 new units\T1\';
% dataDir{2}='W:\Analysis Computer\Persistent Units\Mawake54 new units\T2\';
% dataDir{3}='W:\Analysis Computer\Persistent Units\Mawake54 new units\T3\';
% dataDir{4}='W:\Analysis Computer\Persistent Units\Mawake54 new units\T4\';
dataDir{1}='W:\Analysis Computer\Persistent Units\Mawake64\T1\';
dataDir{2}='W:\Analysis Computer\Persistent Units\Mawake64\T2\';
dataDir{3}='W:\Analysis Computer\Persistent Units\Mawake64\T3\';
dataDir{4}='W:\Analysis Computer\Persistent Units\Mawake64\T4\';

a=[];
b=[];
c=[];
d=[];
base35=[];
base48=[];
for i=1:length(dataDir)
    [curr_a,curr_b,curr_c,curr_d,curr_base35,curr_base48]=calcPersistenceIndex(dataDir{i});
    a=[a; curr_a];
    b=[b; curr_b];
    c=[c; curr_c];
    d=[d; curr_d];
    base35=[base35; curr_base35];
    base48=[base48; curr_base48];
end

figure(); scatter(c./a, d./b);
useUnits=(a>0) & (b>0) & ((c./a)>=0);
disp('total # units');
disp(length(a));
disp('Number of units with a and b not equal to zero');
disp(sum(useUnits));
pers_index=((d./b)-(c./a))./((d./b)+(c./a));
[y,x]=hist(pers_index(useUnits),20);
figure(); plot(x,y);
% return_pers_index=pers_index(useUnits);
return_pers_index=pers_index;

% c2=0.075;
% figure(); scatter(base48(useUnits)-base35(useUnits)+c2.*rand(size(base35(useUnits))),pers_index(useUnits)+c2.*rand(size(base35(useUnits))));
% co=corrcoef(base48(useUnits)-base35(useUnits),pers_index(useUnits));
% disp(co(1,2));

% c2=0.075;
% figure(); scatter(b(useUnits)./a(useUnits)+c2.*rand(size(b(useUnits))),pers_index(useUnits)+c2.*rand(size(b(useUnits))));
% co=corrcoef(b(useUnits)./a(useUnits),pers_index(useUnits));
% disp(co(1,2));

% c2=0.075;
% b_over_a=b(useUnits)./a(useUnits);
% gain_persIndex=(d(useUnits)./b(useUnits))./(c(useUnits)./a(useUnits));
% useuseunits=~isnan(b_over_a) & ~isnan(gain_persIndex) & ~isinf(b_over_a) & ~isinf(gain_persIndex);
% figure(); scatter(b_over_a(useuseunits),gain_persIndex(useuseunits));
% co=corrcoef(b_over_a(useuseunits),gain_persIndex(useuseunits));
% disp(co(1,2));

c2=0.075;
% b_over_a=base48(useUnits)./base35(useUnits);
base_over_a=base48(useUnits)./base35(useUnits);
ev=a(useUnits);
b_over_a=b(useUnits)./a(useUnits);
gain_persIndex=(d(useUnits)./b(useUnits))./(c(useUnits)./a(useUnits));
useuseunits=~isnan(b_over_a) & ~isnan(gain_persIndex) & ~isinf(b_over_a) & ~isinf(gain_persIndex);
figure(); scatter(b_over_a(useuseunits),gain_persIndex(useuseunits));
[co,p]=corrcoef(b_over_a(useuseunits),gain_persIndex(useuseunits));
disp(co(1,2));
disp(p);

sub_gain_persIndex=gain_persIndex(useuseunits);
sub_ev=ev(useuseunits);
sub_b_over_a=b_over_a(useuseunits);
figure();
cmap=colormap(jet);
colmap=cmap(1:8:end,:);
steps=[min(sub_gain_persIndex)-1 0 1 2 3 4 6 8 max(sub_gain_persIndex)+1];
for i=1:length(sub_gain_persIndex)
    sind=find(steps>sub_gain_persIndex(i),1,'first');
    scatter(sub_ev(i),sub_b_over_a(i),[],colmap(sind-1,:));
    hold on;
end