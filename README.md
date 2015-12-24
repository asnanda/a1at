# a1at
functions for dealing with CW EPR spectra of a1at mutants

Changes
========

24/12/15
========

Added examples for generating delH plots, cleaned up some code, and struct2mat_mutant is no longer needed.


10/12/15
=========
You can find examples for making a panel in fourbtwelve_figure.m or Example_script.m

The generate_struct function can now take muultiple directores as inputs,
Or open a UI window for you to select them.

> Ex. files = generate_struct('/Users/Data1','/Users/Data2');


The mutant_struct2mat function has been renamed struct2mat_mutant to keep it consistent

> Ex. mutantnames = unique(struct2mat_mutant(files,'mutant'));


plot_mutant now has some new arguements that can help with cleaning up plots
+ The file struct
+ The residues to be graphed, defaults to 'all'
+ The style of figure. Right now two options, 'default' which look like our common EPR figures, or 'figure' which looks like the panel
+ The color scheme selected. I've preloaded a few in the plot_mutant file, add your own as requierd. Current ones are 
..'default',
..'4b12figure'
..'5e3figure', 
..'dualbinding'
.. Alternatively you can specificy your own colors for specific spectra types (i.e monomerApo)

> Ex. plot_mutant(files,'32','figure','4b12figure');
Here, the configuration of 'residue', 'figure', and '4b12figure' are used.


An alternative might be:
> Ex. plot_mutant(files,'32','figure',{'monomerApo','g','monomer4b12','m'});


