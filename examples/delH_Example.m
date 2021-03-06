clc
clear all
close all

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');

fig = figure(1);
%title('SSM^{-1} vs  \DeltaH^{-1} for Monomer,4b12')

x1 = [files.delH];
x2 = [files2.delH];
y1 = [files.ssm];
y2 = [files2.ssm]; 

hold on
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
scatter(x1,y1,'k','filled');
scatter(x2,y2,'r','filled');
axis([1 8 -12 -7.5])

ax = fig.CurrentAxes.Position;
ax_lims = [fig.CurrentAxes.XLim fig.CurrentAxes.YLim];

magnitude = [];

for m = 1:length(files)
    for z = 1:length(files2)
        if strcmp(files(m).mutant, files2(z).mutant)
            coords = [x1(m) y1(m) x2(z) y2(z)]; %normalize all the coordiantes so that they can be graphed with arrow function
            magnitude(end+1) = sqrt((coords(3)-coords(1)).^2+(coords(4)-coords(2)).^2);
            coords(1) = ax(1)+ ((coords(1)-ax_lims(1))/(ax_lims(2)-ax_lims(1)))*ax(3);
            coords(2) = ax(2)+ ((coords(2)-ax_lims(3))/(ax_lims(4)-ax_lims(3)))*ax(4);
            coords(3) = ax(1)+ ((coords(3)-ax_lims(1))/(ax_lims(2)-ax_lims(1)))*ax(3);
            coords(4) = ax(2)+ ((coords(4)-ax_lims(3))/(ax_lims(4)-ax_lims(3)))*ax(4);
            annotation('arrow',[coords(1) coords(3)],[coords(2) coords(4)],'HeadStyle','vback1'); % only accepts normalized corodiantes
        end
    end
end


hold off
legend('Monomer','4b12');
xlabel('\DeltaH^{-1}','fontsize',14)
ylabel('log(SSM)','fontsize',14);

labelpoints(x1,y1,[files.mutant],'NW','Fontsize',12) %sufficent to only label one with arrows.









