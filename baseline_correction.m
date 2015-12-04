clc
clear all
close all


files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/all/');


figure(1)
mutantnames = unique(mutant_struct2mat(files,'mutant'));
for s = 1:12
    subplot(4,3,s)
    plot_mutant(files,num2str(mutantnames(s)),0);
end
%legend('polymer','pol5e3')



break

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/all/');
%files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Polymer/');
%files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Pol5e3/');

%second spectral moment^{-1} vs
figure(1)
title('\DeltaH^{-1} vs  Residue / Monomer, Monomer4b12')
set (gcf,'PaperUnits','centimeters')
set (gcf,'Position',[0,0,1000,800])
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
set(gca,'xticklabel',(sort(mutant_struct2mat(files,'mutant'))))
set(gca,'xtick',(sort(mutant_struct2mat(files,'mutant'))),'FontName','Arial')
grid off

hold on
bar(mutant_struct2mat(files,'delH'),mutant_struct2mat(files,'ssm'),'k')
scatter((mutant_struct2mat(files2,'delH')),mutant_struct2mat(files2,'ssm'),'r','filled');
scatter((mutant_struct2mat(files3,'delH')),mutant_struct2mat(files3,'ssm'),'b','filled');

hold off
legend('Polymer','Pol5E3');
xlabel('Residue Number','fontsize',12)
ylabel('\DeltaH^{-1}','fontsize',12);




%labelpoints(mutant_struct2mat(files,'delH'),mutant_struct2mat(files,'ssm'),mutant_struct2mat(files,'mutant'),'N')
labelpoints(mutant_struct2mat(files2,'delH'),mutant_struct2mat(files2,'ssm'),mutant_struct2mat(files2,'mutant'),'E')
labelpoints(mutant_struct2mat(files3,'delH'),mutant_struct2mat(files3,'ssm'),mutant_struct2mat(files3,'mutant'),'S')


break

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
%files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5E3DATA/Polymer/'); 
%files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5E3DATA/Polymer5E3/'); 
%files4 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5E3DATA/Monomer5E3/'); 

hold on
scatter(mutant_struct2mat(files,'delH'),mutant_struct2mat(files,'ssmv2'),'k');
%scatter(mutant_struct2mat(files2,'delH'),mutant_struct2mat(files2,'ssmv2'),'r');
%scatter(mutant_struct2mat(files3,'mutant'),mutant_struct2mat(files3,'delH'),'b');



labelpoints(mutant_struct2mat(files,'delH'),mutant_struct2mat(files,'ssm'),mutant_struct2mat(files,'mutant'),'NE')
%labelpoints(mutant_struct2mat(files2,'delH'),mutant_struct2mat(files2,'ssm'),mutant_struct2mat(files2,'mutant'),'SE')
%labelpoints(mutant_struct2mat(files3,'mutant'),mutant_struct2mat(files3,'delH'),mutant_struct2mat(files3,'antibody'),'E')


