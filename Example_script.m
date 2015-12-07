clc
clear all
close all

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');

fig = figure(1);
title('\DeltaH^{-1} vs  \DeltaH^{-1}  Monomer,4b12')

x1 = (struct2mat_mutant(files,'delH'));
x2 = (struct2mat_mutant(files2,'delH'));
y1 = struct2mat_mutant(files,'ssm'); 
y2 = struct2mat_mutant(files2,'ssm');

hold on
scatter(x1,y1,'k','filled');
scatter(x2,y2,'r','filled');
ax = fig.CurrentAxes.Position;
ax_lims = [fig.CurrentAxes.XLim fig.CurrentAxes.YLim];


% arrow generates -> only works for two, but coudl be generalized. 
for m = 1:length(files)
    for z = 1:length(files2)
        if strcmp(files(m).mutant, files2(z).mutant)
            coords = [x1(m) y1(m) x2(z) y2(z)]; %normalize all the coordiantes so that they can be graphed with arrow function
            coords(1) = ax(1)+ ((coords(1)-ax_lims(1))/(ax_lims(2)-ax_lims(1)))*ax(3);
            coords(2) = ax(2)+ ((coords(2)-ax_lims(3))/(ax_lims(4)-ax_lims(3)))*ax(4);
            coords(3) = ax(1)+ ((coords(3)-ax_lims(1))/(ax_lims(2)-ax_lims(1)))*ax(3);
            coords(4) = ax(2)+ ((coords(4)-ax_lims(3))/(ax_lims(4)-ax_lims(3)))*ax(4);
            annotation('arrow',[coords(1) coords(3)],[coords(2) coords(4)]) % only accepts normalized corodiantes
        end
    end
end



hold off
legend('Monomer','4b12');
xlabel('\DeltaH^{-1}','fontsize',12)
ylabel('SSM^{-1}','fontsize',12);

labelpoints(x1,y1,struct2mat_mutant(files,'mutant'),'NW') %sufficent to only label one with arrows.
%labelpoints(x2,y2,struct2mat_mutant(files2,'mutant'),'S')
%labelpoints(mutant_struct2mat(files3,'delH'),mutant_struct2mat(files3,'ssm'),mutant_struct2mat(files3,'mutant'),'S')










