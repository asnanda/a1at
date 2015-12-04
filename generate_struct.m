
function files = generate_struct(varargin)
%define the empty struct -> files
files = struct('filepath',{},'name',{},'mutant',{},'state',{},'antibody',{},'x',{},'y',{},'si',{},'di',{},'di_val',{},'y_cor',{},'si_cor',{},'di_val_cor',{},'spectra_params',{},'y_cor_norm',{},'delH',{},'fsm',{},'ssm',{});

% either load files, or open them.
switch nargin
    case 0
        % select the files
        [filenames,filepath] = uigetfile('*.spc','Bruker SPC Files Only (*.spc)','MultiSelect','on');
        filenames = unique(filenames);
    case 1
        % handle a directory
        filepath = varargin{1};
        tmp_filenames = dir(strcat(filepath,'*.spc'));
        filenames = mutant_struct2mat(tmp_filenames,'name');
        
    otherwise
        warning('Cannot handle more than one file');
        % filenames must have the spc on it. shoud be filepath, filename,
        % sc
        
        
    
    otherwise
        
        
        
        warning('Not ready to handle manual arguements! Ahhh!')
        return
       
        
        % how to handle mutiple file_paths
            
         % check all files are unique 
        
        
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
     % baseline corretion is carried out by the bc_polyfit function, found
     % in bc_polyfit.m
     % current paramters used are "linear" for EPR correction, and "Poly2"
     % for absorbance spectra correction. 
     
     files(i).y_cor = bc_polyfit(files(i).x,files(i).y,'epr','linear');
     files(i).si_cor = bc_polyfit(files(i).x,cumsum(files(i).y_cor),'absorbance','poly2');
     files(i).di_val_cor = trapz(files(i).x,files(i).si_cor);
     files(i).y_cor_norm = files(i).y_cor/files(i).di_val_cor;
    
    % generate spectral summary measurements--------------------------------
    files(i).delH = 1/(files(i).x(files(i).y == min(files(i).y)) - files(i).x(files(i).y == max(files(i).y)));
    files(i).fsm = moment(files(i).y_cor_norm,1);
    files(i).fsmv2 = 1/get_mutant_moment(files,i,1);
    files(i).ssm = 1/(moment(files(i).y_cor_norm,2));
    files(i).ssmv2 = 1/get_mutant_moment(files,i,2);
     
     
    

end % loop through files end





end % generate_struct end





function result = get_characteristic(characteristic,name)
switch characteristic
    case 'name'
        [~,name,~] = fileparts(char(name));
        result = name;
        
    case 'state'
        if length(regexp(name,'pol','match')) >= 1
            result = {'polymer'}; 
        else
            result = {'monomer'};
        end

    case 'antibody'
        if length(regexpi(name,'4b12','match')) >= 1
            result = {'4b12'};
        elseif length(regexpi(name,'5e3','match')) >= 1
            result = {'5e3'};
        else
            result = {'Apo'};
        end
        
    case 'mutant'
        match = regexp(name,'[A-Z](\d+)[A-Z]|(\d+)[C]|(\d+)','match');
        if length(match) >= 1
            match2 = regexp(match(1),'(\d+)','match'); % match(1) is used because the resisude number comes before any other numbers in a file name
            result = match2{1};
        else
            result = name; 
        end

end % end switch statement
end % end funciton


function result = get_mutant_moment(files,n,order)

switch order
    case 1
        x_vals = files(n).x;
        y_vals = x_vals'.*files(n).y_cor_norm;
        fsm = trapz(y_vals)/trapz(files(n).y_cor_norm);
        result = fsm;
       
    case 2
        fsm = get_mutant_moment(files,n,1);
        x_vals = (files(n).x - fsm).^2;
        y_vals = x_vals'.*files(n).y_cor_norm;
        ssm = trapz(y_vals)/trapz(files(n).y_cor_norm);
        result = ssm;
end % end N switch
end



% generate_struct: Open .spc files from our CW EPR experiments and
% load the data, and characteristics, into a struct for future manipluation
%
% generate_struct: when run without any inputs, opens a GUI so that the we can
% select epr files to open. It can also accept a path to a file as
% an input if the path is put in 'quotes' with the extension '.spc'
%
% generate_struct -> UI. window
% generate_struct('path1','path2','path3'...)
%
% struct files = generate_struct
% struct files = generate_struct('path1','path2','path3'...) 


% generate_struct: A function that opens spc files related to the antibody 
%CW experiments and extracts their data for further mainpulation.
%Information pulled includes: filename, file path, antibody, state,
%magnetic field data, intensity data, uncorrected single integral values,
%uuncorrected double integral values, uncorreted singular value double
%integral,corrected single integral values,
%corrected double integral values, corrected singular value double
%integral,
%
%
% Inputs:
%    input0     - none, results in a UI for selection
%    inputn     - string input to the path of a file
%
% Outputs:
%    output0    - a struct, files, storing the following data:

% Example: 
%    files = generate_struct
%    files = generate_struct('/path/to/file.spc', /another/path/to/files.scp,....')
%          
%
% Other m-files required:   bc_polyfit.m - for polynomial basline
%                           correction
%
% Subfunctions:             get_characteristics
%
% MAT-files required:       easyspin package for eprload
%


