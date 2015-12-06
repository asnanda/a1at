clc
clear all
close all

%CHANGES

% the generate_struct function can now take muultiple directores as inputs,
% or open a UI window for you to select them.
files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/','/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');

mutantnames = unique(struct2mat_mutant(files,'mutant'));

figure(1)

for s = 1:6
    subplot(3,2,s)
    plot_mutant(files,num2str(mutantnames(s)),'figure','4b12figure');
end

xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
