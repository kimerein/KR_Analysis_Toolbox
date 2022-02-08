function copyFig1ToFig2()

first_fig = figure(1);
second_fig = figure(2);
first_ax = findobj(first_fig, 'type', 'axes');
second_ax = findobj(second_fig, 'type', 'axes');
if length(first_ax) ~= 1 || length(second_ax) ~= 1
error('this code requires the two figures to have exactly one axes each');
end
ch2 = get(second_ax, 'children'); %direct children only and don't try to find the hidden ones
copyobj(ch2, first_ax); %beam them over

end