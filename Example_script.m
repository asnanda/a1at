clc
clear all
close all

files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
%second spectral moment^{-1} vs

figure(1)
title('\DeltaH^{-1} vs  Residue / Monomer, Monomer4b12')

hold on
scatter((struct2mat_mutant(files,'delH')),struct2mat_mutant(files,'ssmv2'),'k','filled');
scatter((struct2mat_mutant(files2,'delH')),struct2mat_mutant(files2,'ssmv2'),'r','filled');


for m = 1:length(files)
    for z = 1:length(files2)
        if strcmp(files(m).mutant, files2(z).mutant)
            quiver(files(m).delH,files(m).ssmv2,files2(z).delH -files(m).delH,files2(z).ssmv2 - files(m).ssmv2)
        end
    end
end



hold off
legend('Monomer','4b12');
xlabel('\DeltaH^{-1}','fontsize',12)
ylabel('SSM^{-1}','fontsize',12);

labelpoints(struct2mat_mutant(files,'delH'),struct2mat_mutant(files,'ssmv2'),struct2mat_mutant(files,'mutant'),'N')
labelpoints(struct2mat_mutant(files2,'delH'),struct2mat_mutant(files2,'ssmv2'),struct2mat_mutant(files2,'mutant'),'S')
%labelpoints(mutant_struct2mat(files3,'delH'),mutant_struct2mat(files3,'ssm'),mutant_struct2mat(files3,'mutant'),'S')





break


%CHANGES

% the generate_struct function can now take muultiple directores as inputs,
% or open a UI window for you to select them.
files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/all/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Polymer/');
files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Monomer/','/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Polymer/');

% the mutant_struct2mat function has been renamed struct2mat_mutant to keep
% it consistent
mutantnames = unique(struct2mat_mutant(files,'mutant'));

figure(1)

for s = 1:9
    subplot(3,3,s)
    %plot mutant now has some new arguements that can help
    % 1) The file struct
    % 2) The residues to be graphed, defaults to 'all'
    % 3) The style of figure. Right now two options, 'default' which looks
    % like our common EPR figures, or 'figure' which looks like the panel
    % 4) The color scheme selected. I've preloaded a few in the plot_mutant
    % file, add your own as requierd. Current ones are 'default',
    % '4b12figure', '5e3figure', 'dualbinding'. Alternatively you can
    % specificy your own colors for specific spectra types (i.e monomer
    % Apo)
    % Here, the configuration of 'residue', 'figure', and '4b12figure' are
    % used to generate a panel for the 4b12 mutants
    
    plot_mutant(files,num2str(mutantnames(s)),'figure','4b12figure');
    
    % an alternative might be:
    % plot_mutant(files,num2str(mutantnames(s)),'figure',{'monomerApo','g','monomer4b12','m'});
end







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


