function  result = struct2mat_mutant(S,field)

switch ischar(S(1).(deblank(field))) % if the field is a character array, then output a string.
    case 0 % else
        switch iscell(S(1).(deblank(field))) % check if it's a cell-string
            case 0 % if it isnt, then output the numerical value
                 for n = 1:length({S.(deblank(field))})
                    struct_sum(n) = S(n).(deblank(field));
                 end
            case 1 % if it is, then convert as required
                switch field
                    case 'antibody'
                        for n = 1:length({S.(deblank(field))})
                            struct_sum(n) = cellstr(S(n).(deblank(field)));
                        end
                    case 'state'
                        for n = 1:length({S.(deblank(field))})
                            struct_sum(n) = cellstr(S(n).(deblank(field)));
                        end
                    otherwise % mutant
                        for n = 1:length({S.(deblank(field))})
                            struct_sum(n) = str2double(S(n).(deblank(field)));
                        end
                end               
        end
       
    case 1
        for n = 1:length({S.(deblank(field))})
            struct_sum(n) = cellstr(S(n).(deblank(field)));
        end
end


result = struct_sum;
end
