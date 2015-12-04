function  result = mutant_struct2mat(S,field)


switch ischar(S(1).(deblank(field)))
    case 0
        switch iscell(S(1).(deblank(field)))
            case 0
                 for n = 1:length({S.(deblank(field))})
                    struct_sum(n) = S(n).(deblank(field));
                 end
            case 1
                switch field
                    case 'antibody'
                        for n = 1:length({S.(deblank(field))})
                            struct_sum(n) = cellstr(S(n).(deblank(field)));
                        end
                    case 'state'
                        for n = 1:length({S.(deblank(field))})
                            struct_sum(n) = cellstr(S(n).(deblank(field)));
                        end
                    otherwise
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
