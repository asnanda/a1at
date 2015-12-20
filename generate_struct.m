% generate_struct: Open .spc files from our CW EPR experiments and
% load the data, and characteristics, into a struct for  manipulation.
%
% generate_struct: when run without any inputs, opens a GUI so that the you can
% select epr files to open. It can also accept a path to a directory
% as an input if the path is put in 'quotes'.
%IN
%        generate_struct;  -> UI. window
%        generate_struct('directory1')
%        generate_struct('directory1','directory2'..etc) NOTE -> Inefficent
%OUT
%        struct files = generate_struct;
%        struct files = generate_struct('path1','path2','path3'...); 
%
%
% generate_struct: A function that opens spc files related to the antibody 
% CW experiments and extracts their data for further mainpulation. It pulls
% relevant spectral parameters, generates normalized and baseline
% corrected intensity measurements, and generates calculated spectral
% values (delH, spectral moments, etc)
%
%
% Inputs:
%    input0     - none, results in a UI for selection
%    input1     - string input to the path of a file
%
% Outputs:
%  output0  - a struct  storing the following data:
%     file name                        .name             ex. {'name'}
%     file path                        .filepath         ex. {'file/path/file.spc'}
%     mutant                           .mutant           ex. {'32'}
%     antibody                         .antibody         ex. {'Apo'} / {'4b12'} / {'5e3'} 
%     state                            .state            ex. {'monomer'} / {'polymer'}
%     magnetic field data (mT)         .x                ex. [val1...valn]
%     intensity data                   .y                ex. [val1...valn]
%     uncor. single integral           .si               ex. [val1....valn]
%     uncor. double integral value     .di_val           ex. (2e8)
%     cor. intensity data,             .y_cor            ex. [val1...valn]
%     cor. singular integral           .si_cor           ex. [val1...valn]
%     cor. double integral value       .di_val_cor       ex. (2e8)
%     normalized cor. intensity data   .y_cor_norm       ex. [val1....valn]
%     spectral parameters              .spectra_params   ex. [MWFreq,Npoints..etc]
%     1/?H,center peak separation(mT)  .delH             ex. (4.223)
%     1/first CENTRAL spec. moment     .fsm              ex. (0)
%     1/second CENTRAL spec. moment    .ssm              ex. (4.0e33)

% Example: 
%    files = generate_struct;
%    files = generate_struct('/directory/with/files/')
%          
%
% Other m-files required:   baseline_correct_mutant.m - for  basline correction
%
% Subfunctions:             get_characteristics(characteristic, name)
%                               - parses file names,
%                               - pulls trial conditions (antibody, state, etc)
%                           mutant_struct_sort(struct)
%                               - sorts struct by first by mutant, then
%                               state, then antibody. Mutants by increasing number,
%                               monomer before polymer, Apo before 4b12/5e3, and  
%                           get_mutant_moment
%                               - attempted implementation of the first
%                               spectral moment, current not working
%
% MAT-files required:       struct2mat_mutant.m
%                           basline_correct_mutant.m
%                           easyspin package for eprload
%

function files = generate_struct(varargin)
%define the empty struct -> files
files = struct('filepath',{},'name',{},'mutant',{},'state',{},'antibody',{},'x',{},'y',{},'si',{},'di',{},'di_val',{},'y_cor',{},'si_cor',{},'di_val_cor',{},'y_cor_norm',{},'spectra_params',{},'delH',{},'fsm',{},'ssm',{});

% either load files, or open them.
switch nargin
    case 0
        % select the files
        [filenames,filepath] = uigetfile('*.spc','Bruker SPC Files Only (*.spc)','MultiSelect','on');
        filenames = unique(filenames);
    case 1
        % handle a directory
        filepath = varargin{1}; % get the first arguement
        tmp_filenames = dir(strcat(filepath,'*.spc')); % get the spc files
        filenames = struct2mat_mutant(tmp_filenames,'name'); % extract filenames from the 'dir' struct. 
        
    otherwise % more than one directory (uses recursion) -> can take time
        % call the generate struct function on the directories passed in,
        first = generate_struct(varargin{1});
        for r = 2:nargin
            first = cat(2,first,generate_struct(varargin{r})); % join them together into one struct
        end
        
        files = mutant_struct_sort(first);
       
        return
        
        
end % nargin end

%--------------------------------------------------------------------------

% define variables

NUMBER_OF_FILES = length(filenames); % get the number of files selected

% load the struct
for i = 1:NUMBER_OF_FILES
    % get file ------------------------------------------------------------
    files(i).filepath = strcat(char(filepath),char(filenames(i))); 
    % get name ------------------------------------------------------------
    files(i).name = get_characteristic('name',filenames(i));
    % get mutant ----------------------------------------------------------
    files(i).mutant = get_characteristic('mutant',files(i).name);
    % get state ----------------------
    files(i).state = get_characteristic('state',files(i).name);
    % get antibody --------------------------------------------------------
    files(i).antibody = get_characteristic('antibody',files(i).name);
    % get raw/uncorrected magnetic field, intensity data-------------------
    try
        [b,s,Params] = eprload(char(files(i).filepath)); % eprload from easyspin
    catch
        warning('No files loaded');
        break 
    end

    b = b/10; %convert to millitesla
    files(i).x = b'; % fix alignment error, make ararys column vectors by transpose
    files(i).y = s;
    files(i).si = cumsum(s); % cumsum to generate the single integral
    files(i).di = cumsum(cumsum(s)); 
    files(i).di_val = trapz(cumsum(s)); % a value for the area under the absorbance curve
    files(i).spectra_params = Params;
    
    % get corrected magnetic field, intensity data------------------------
    % baseline corretion is carried out by the baseline_correct_mutant function
    % found in baseline_correct_mutant.m
    % current paramters used are "linear" for EPR correction, and "poly2"
    % for absorbance spectra correction. 
     
    files(i).y_cor = baseline_correct_mutant(files(i).x,files(i).y,'epr','linear');
    files(i).si_cor = baseline_correct_mutant(files(i).x,cumsum(files(i).y_cor),'absorbance','poly2');
    files(i).di_val_cor = trapz(files(i).x,files(i).si_cor);
    files(i).y_cor_norm = files(i).y_cor/files(i).di_val_cor;
    files(i).y_cor_formoment = files(i).y_cor/trapz(files(i).y_cor);
    files(i).si_cor_norm = files(i).si_cor / max(files(i).si_cor);
    files(i).si_cor_norm_di = files(i).si_cor/ (files(i).di_val_cor);
    files(i).si_cor_norm_si = files(i).si_cor / trapz(files(i).si_cor);
        
    % generate spectral summary measurements--------------------------------
    files(i).delH = 1/(files(i).x(files(i).y == min(files(i).y)) - files(i).x(files(i).y == max(files(i).y)));
    files(i).p2pA = max(files(i).y_cor_norm)-min(files(i).y_cor_norm);
    
    files(i).fsm = moment(files(i).y_cor_norm,1); % should always be zero
    files(i).ssm = log(moment(files(i).y_cor_norm,2));

    files(i).fsmv2 =get_mutant_moment(files,i,1);
    files(i).ssmv2 = log(1/get_mutant_moment(files,i,2));

end % loop through files end



% sort the struct in the following way:
% Mutant(numerical), state(mon, pol), antibody(apo, 4b12, 5e3)

files = mutant_struct_sort(files);



end % generate_struct end



% main function for parsing file names
function result = get_characteristic(characteristic,name)
switch characteristic
    case 'name'
        [~,name,~] = fileparts(char(name));
        result = name;
        
    case 'state' % parse for state
        if length(regexp(name,'pol','match')) >= 1 % if match for 'pol'
            result = {'polymer'}; % -> set as polymer
        else
            result = {'monomer'}; %-> else set as monomer
        end

    case 'antibody'
        if length(regexpi(name,'4b12','match')) >= 1 % if a hit for 4b12,
            if length(regexpi(name,'5e3','match')) >= 1 % check for dual 5e3, handles dual binding condition
                result = {'4b12-5e3'};
            else 
                result = {'4b12'}; %else set as 4b12
            end
           
        elseif length(regexpi(name,'5e3','match')) >= 1 %check for hit for 5e3
            result = {'5e3'};
        else
            result = {'Apo'}; % if no hits, then -> Apo
        end
        
    case 'mutant'
        match = regexp(name,'[A-Z](\d+)[A-Z]|(\d+)[C]|(\d+)','match'); %regex for matching '123C' type names
        if length(match) >= 1
            match2 = regexp(match(1),'(\d+)','match'); % match(1) is used because the resisude number comes before any other numbers in a file name
            result = match2{1};
        else
            result = name; 
        end

end % end switch statement
end % end function


%function for sorting struct array according to "antibody_order" variable
function files_sorted = mutant_struct_sort(files)
disp('Loaded');
mutes = unique(struct2mat_mutant(files,'mutant'));
files_sorted = files(1);
antibody_order = {'monomerApo','monomer4b12','monomer5e3','monomer4b12-5e3','monomer5e3-4b12','polymerApo','polymer4b12','polymer5e3','polymer4b12-5e3','polymer5e3-4b12'};
hits = [];
for i = 1:length(mutes) % walk through the list of unique mutants
    for j = 1:length(files) % walk through the struct.
        if str2double(files(j).mutant) == mutes(i) % if the mutants line up,
            hits(end+1) = j; % store the indexes of those mutant sin 'hits'
        end
    end
    
    for k = 1:length(antibody_order) % walk through the order 
        for l = 1:length(hits) % walk through 'hits' generated above for a specific mutant.
            if strcmp(strcat(files(hits(l)).state,files(hits(l)).antibody),antibody_order(k)) % if the entry in the order matches the entry in files,
                files_sorted(end+1) = files(hits(l)); % then add it to the 'files_sorted' struct.
            end
        end
    end
    hits = []; %clean the hits variable, so it can be reused for the next mutant.
end
files_sorted = files_sorted(2:length(files_sorted)); %delete the value used to intialize the array
return
end


% new implementation to calcualte ssm, fsm. Not working.
function result = get_mutant_moment(files,n,order)

switch order
    case 1
        result = 0; % doesnt' involve calcuating the 1st moment
    case 2
        result = (trapz(((files(n).x).^2)'.*files(n).y_cor_norm)/trapz(files(n).y_cor_norm)) - (trapz((files(n).x)'.*files(n).y_cor_norm)/trapz(files(n).y_cor_norm)).^2;
end


return
%{
%first attempt
switch order
    case 1
        x_vals = files(n).x;
        y_vals = x_vals'.*files(n).y_cor_formoment;
        fsm = trapz(y_vals);
        result = fsm;
       
    case 2
        fsm = get_mutant_moment(files,n,1);
        x_vals = (files(n).x - fsm).^2;
        y_vals = x_vals'.*files(n).y_cor_formoment;
        ssm = trapz(y_vals);
        result = ssm;
end % end N switch
%}
end




