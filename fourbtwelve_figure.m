clc
clear all
close all
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
