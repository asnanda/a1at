% working. 






function [corrected_value1, corrected_value2] = baseline_correct_mutant(b,s,spectra,fit_type)
if isempty(s) == 0
    switch spectra
        case 'epr'
            DATA_SET_SIZE = 50;  % large window, due to larger number of points
            [up_start, up_end] = get_indicies(s,'epr','upfield',DATA_SET_SIZE);
            [down_start, down_end] = get_indicies(s,'epr','downfield',DATA_SET_SIZE);
            
            y_set1 = cat(1,(s(up_start:up_end)), s(down_start:down_end)); % generated with multiple data points
            x_set1 = cat(1,(b(up_start:up_end)), b(down_start:down_end));
                    
            y_set2 = [median(s(up_start:up_end)), median(s(down_start:down_end))]; %generated w medians
            %prevents overlaps.
            x_set2_1 = b(s == y_set2(1));
            x_set2_2 = b(s == y_set2(2));
            x_set2 = [x_set2_1(1), x_set2_2(1)]; 
            
            switch fit_type
                case 'poly'
                    
                    [poly_coeff1,Stc] = polyfit(x_set1,y_set1,2);
                    %[poly_coeff2,Stc] = polyfit(x_set2,y_set2,2); badly conditioned, needs a third arguement

                   
                    for n = 1:length(s)
                        corrected_value1(n) = s(n) - (poly_coeff1(1)*(b(n)^2)+poly_coeff1(2)*b(n)+poly_coeff1(3));
                        corrected_value2(n) = 0; %s(n) - (poly_coeff2(1)*(b(n)^2)+poly_coeff2(2)*b(n)+poly_coeff2(3));
                        % not returned, because badly conditioned with only
                        % two points. instead, multiple point set fit
                        % polynomial is returned. Not reccomended.
                    end
                    
                otherwise % linear, catch all
                    
                    [lin_coeff1,Stc] = polyfit(x_set1,y_set1,1); % now with MDP
                    [lin_coeff2,Stc] = polyfit(x_set2,y_set2,1); % now the median
                    

                    for n = 1:length(s)
                        corrected_value1(n) = s(n) - (lin_coeff1(1)*(b(n))+lin_coeff1(2)); % MDP 
                        corrected_value2(n) = s(n) - (lin_coeff2(1)*(b(n))+lin_coeff2(2)); % Median
                    end        
                    
            end %fit_type switch end
            
        
        case 'absorbance'
            %smaller window due to smaller corrections -> polyfit. 
            % note here, because of the applicaiton of cumsum to the orignal values
            % in order to achieve the absorbance spectra from teh corrected
            % y values:
            % the dimensions of the arrays change. its easier to flip the
            % y_set arrays back to rows x 1 column. 
            % This transform only need be applied to the MDP sets, the data
            % extracted directly from the array, rather than being wrapped
            % around it.
            
            
            DATA_SET_SIZE = 5;
            [up_start, up_end] = get_indicies(s,'absorbance','upfield',DATA_SET_SIZE);
            [down_start, down_end] = get_indicies(s,'absorbance','downfield',DATA_SET_SIZE);
            
            y_set1 = cat(2,(s(up_start:up_end)), s(down_start:down_end));% generated with multiple data points
            x_set1 = cat(1,(b(up_start:up_end)), b(down_start:down_end));
         
                    
      
           
            switch fit_type
                case 'linear'
                    warning('Warning: Yields a poor fit, and a poor correction.')
                    [lin_coeff1,Stc] = polyfit(x_set1,y_set1',1);
                    
                    for n = 1:length(s)
                        corrected_value1(n) = s(n) - (lin_coeff1(1)*(b(n))+lin_coeff1(2));
                        corrected_value2(n) = 0; % no secondary value is calculated.
                    end
                    
                
                case 'poly2'
                    % empircal evidence suggests that downfield paramters need to
                    % be between the lowest -> 344 mT to get best fit results.
                    % averaging can continue at the high field, but with a higher
                    % paramter window, maybe 20
                    
                     y_set1 = cat(2,s(end), s(1:75));% generated with multiple data points
                     x_set1 = cat(1,b(end), b(1:75));
         
                    [poly_coeff1,Stc] = polyfit(x_set1,y_set1',2);

                    for n = 1:length(s)
                        corrected_value1(n) = s(n) - (poly_coeff1(1)*(b(n)^2)+poly_coeff1(2)*b(n)+poly_coeff1(3));
                        corrected_value2(n) = 0; % no secondary value is put out.
                    end
                    
                    
                
                otherwise % normally returns poly 
                    
                    [poly_coeff1,Stc] = polyfit(x_set1,y_set1',2);

                    for n = 1:length(s)
                        corrected_value1(n) = s(n) - (poly_coeff1(1)*(b(n)^2)+poly_coeff1(2)*b(n)+poly_coeff1(3));
                        corrected_value2(n) = 0;
                    end

            end %fit_type switch end
                    

        otherwise
            warning('spectra type not indicated, either epr, or absorbance. Exiting')
            return
    end 

 % check if the point set is empty, if it is, skip to else statement and exit.
    
% generate downfield point set
% select multipe point sets, and find the median
% if the values 'agree', then set point set as downfied. 
% sample sets that minimze variaton .



%generate upfield point set

% do these generaitons via asemlbing a distribution of the values, finding
% the mean 


else
    warning('Error calculating baseline, point set returned');
    corrected_value1 = s;
    
end % end outer protective loop


end



function [start_index,end_index] = get_indicies(s,type,field,DATA_SET_SIZE)
switch type
    case 'epr'
        switch field
            case 'upfield'
                window_start = length(s);
                window_end = length(s)-100+DATA_SET_SIZE;
                step_direction = -1;

            case 'downfield'
                window_start = 2;
                window_end = 102 - DATA_SET_SIZE;
                step_direction = 1;

            otherwise
                window_start = 2;
                window_end = 100 - DATA_SET_SIZE;
        end
        
        pre_variance = 2e40; %define the initial variance parameter 
        for i = window_start:step_direction:window_end
            j = i+step_direction*DATA_SET_SIZE; 
            if i < j
                variance = var(s(i:j));
                 if variance < pre_variance
                    start_index = i;
                    end_index = j;
                    pre_variance = variance;
                end
                
            else
                variance = var(s(j:i));
                if variance < pre_variance
                    start_index = j;
                    end_index = i;
                    pre_variance = variance;
                end
            end
        end
        
    case 'absorbance'
        switch field
            case 'upfield'
                window_start = length(s);
                window_end = length(s)-10+DATA_SET_SIZE;
                step_direction = -1;

            case 'downfield'
                window_start = 2;
                window_end = 12 - DATA_SET_SIZE;
                step_direction = 1;

            otherwise
                window_start = 2;
                window_end = 10 - DATA_SET_SIZE;
        end
        
        pre_variance = 2e40; %define the initial variance parameter 
        for i = window_start:step_direction:window_end
            j = i+step_direction*DATA_SET_SIZE; 
            if i < j
                variance = var(s(i:j));
                 if variance < pre_variance
                    start_index = i;
                    end_index = j;
                    pre_variance = variance;
                end
                
            else
                variance = var(s(j:i));
                if variance < pre_variance
                    start_index = j;
                    end_index = i;
                    pre_variance = variance;
                end
            end
        end

        
        
        
        
end % end switch type
end % end get_indicies function


