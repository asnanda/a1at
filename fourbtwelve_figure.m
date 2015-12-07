clc
clear all
close all

%CHANGES

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/','/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
front = imread('/Users/Scott/Desktop/UCL Lab/Graphics/JournalFigures/front.png');


mutantnames = unique(struct2mat_mutant(files,'mutant'));
resnames =  {'E','S','H','N','E','L','V','T','I','T'}

figure(1)
set (gcf,'PaperUnits','centimeters')
set (gcf,'Position',[0,0,1000,800])
set (gcf,'InvertHardcopy','off','Color',[1 1 1])





offset = 0.08;
totheight = 1-offset;
plot_height = 0.9*(totheight/6);
plot_width = 0.3;


e_front = size(front);
z = 3e-4;

subplot('Position',[0.45 0.5-(z*e_front(1))/2 z*e_front(2)-0.08 z*e_front(1)])
image(front)
set(gca, 'box', 'off','ytick', [],'xtick', [] , 'xcolor', 'w' ,'ycolor', 'w'  )



for i = 1:6
    plot_pos = (6-i)*(totheight/6)+offset;
    subplot('Position',[0.08, plot_pos, plot_width,plot_height])
    plot_mutant(files,num2str(mutantnames(i)),'figure','4b12figure');
    text(0.8,0.6,strcat(resnames(i),num2str(mutantnames(i)),'C') ,'Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');
end
set(gca,'XColor',[0 0 0]);
text(0.3,-.35,'Magnetic Field[mT]','Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');

[ax2,h2]=suplabel('d \chi" / dB','y',0.001);
set(h2,'FontSize',18,'FontName','arial');

