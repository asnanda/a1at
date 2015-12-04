% plot mutant

% spectratype, normalize
function plot_mutant(files,mutant,color_scheme)


if isempty(files)
    warning('No files to plot.')
    return
end

% defaults
%corrected = 'Y';
%mutant = 'all';
%color_scheme = 0;

% set up the color schemes:
switch color_scheme
    case 0
        %        [mon, mon 4b12, pol]
        colors = ['k','b','r'];
    otherwise
        %        [mon, mon 4b12, pol]
        colors = ['k','b','r'];
end % end color scheme switch



switch mutant
    case 'all'
        xzero_standard = get_xcoord_of_zero(files(1).x,files(1).y_cor);
        hold on
        plot_aesthetics;
        for l = 1:length(files)
            x_corr_factor= xzero_standard - get_xcoord_of_zero(files(l).x,files(l).y_cor);
            legend_name = strcat(files(l).mutant,'C-',files(l).state,'-',files(l).antibody);
            plot(files(l).x + x_corr_factor,files(l).y_cor/files(l).di_val_cor,'DisplayName',char(legend_name));
            legend('-DynamicLegend');
            axis([342 355 -inf inf])
            xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
        end
        hold off
        return
        
        
    otherwise
        switch length({mutant})
            otherwise  
                mutant_id = regexp(mutant,'(\d+)','match');
                index_list = [];
                for l = 1:length(files)
                    if strcmp(files(l).mutant,mutant_id) == 1
                        index_list(end+1) = l;
                    end
                end
                if strcmp(mutant_id,'81') == 1
                    a = index_list(1);
                    b = index_list(2);
                    index_list(1) = b;
                    index_list(2) = a;
                end
                
                xzero_standard = get_xcoord_of_zero(files(index_list(1)).x,files(index_list(1)).y_cor);
                hold on
                plot_aesthetics; % call the plot aesthetics.
                for p = 1:length(index_list)
                    cur_index = index_list(p);
                    x_corr_factor= xzero_standard - get_xcoord_of_zero(files(cur_index).x,files(cur_index).y_cor);
                    %legend_name = strcat(files(cur_index).mutant,'C-',files(cur_index).state,'-',files(cur_index).antibdy);
                    legend_name = files(cur_index).antibody;
                    if strcmp(legend_name,'5e3')
                        color = 'r';
                        lspec = '-';
                    else
                        color = 'k';
                        lspec = '-';
                    end
                    plot(files(cur_index).x+x_corr_factor,files(cur_index).y_cor_norm,lspec,'linewidth',1,'color',color,'DisplayName',char(legend_name));
                    legend('-DynamicLegend');
                    axis([343 353 -inf inf])
                    title(strcat(mutant_id,'C'),'FontName','arial','fontweight','normal')
                    xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
                    set(gca,'yticklabel','');
                end
       
                return 
                
                
                
                
                
                % single mutant to graph, align all and graph all
             % for more than one mutant graphed against eachother. 
                %treat each one seperately. 
                
                
        end
        
        
end



end





function xzero = get_xcoord_of_zero(base_x,base_y)
    for t = 1:length(base_y)
        if base_y(t) == max(base_y) % find the maximum
            for r = t:(t+80) % iterate 30 steps forward, usually more than required
                if sign(base_y(r)) == -1 % find where the intensity switches to negative
                    y1 = base_y(r-1);
                    y2 = base_y(r);
                    x1 = base_x(r-1);
                    x2 = base_x(r);
                    mslope = (y2-y1)/(x2-x1);
                    xzero = x1 - (y1/mslope);
                    return
                end
            end
        end 
    end
end


function plot_aesthetics(varargin)
switch nargin
    case 0
        %aesthetics:
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 1024, 800], 'PaperUnits', 'centimeters','PaperSize', [21, 29.7])
        set (gcf,'InvertHardcopy','off','Color',[1 1 1])
        grid off
        box on
        
        
    otherwise
       return
end
end

