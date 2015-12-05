% plot mutant

% spectratype, normalize
function plot_mutant(files,mutant,style,color_scheme)

if isempty(files)
    warning('No files to plot.')
    return
end

if (~exist('mutant','var'))
    mutant = 'all';
end

if (~exist('style','var'))
    style = 'default';
end

if (~exist('color_scheme','var'))
    color_scheme = 'default';
end


if iscell(mutant)
    for z = 1:length(mutant)
        plot_mutant(files,mutant{z},style,color_scheme)
    end
    
else
    switch mutant
        case 'all'
            hold on
            xzero_standard = get_xcoord_of_zero(files(1).x,files(1).y_cor);
            plot_aesthetics(style);
            for l = 1:length(files)
                x_corr_factor= xzero_standard - get_xcoord_of_zero(files(l).x,files(l).y_cor);
                legend_name = strcat(files(l).mutant,'C-',files(l).state,'-',files(l).antibody);
                plot(files(l).x + x_corr_factor,files(l).y_cor/files(l).di_val_cor,'DisplayName',char(legend_name));
                legend('-DynamicLegend');
            end
            hold off
            return
            
        otherwise
            mutant_id = regexp(mutant,'(\d+)','match');
            index_list = [];

            for l = 1:length(files)
                if strcmp(files(l).mutant,mutant_id) == 1
                    index_list(end+1) = l;
                end
            end

            hold on
            xzero_standard = get_xcoord_of_zero(files(index_list(1)).x,files(index_list(1)).y_cor);
            plot_aesthetics(style); % call the plot aesthetics.
            for p = 1:length(index_list)
                cur_index = index_list(p);
                x_corr_factor= xzero_standard - get_xcoord_of_zero(files(cur_index).x,files(cur_index).y_cor);
                %legend_name = strcat(files(cur_index).mutant,'C-',files(cur_index).state,'-',files(cur_index).antibdy);
                legend_name = files(cur_index).antibody;
                [lspec,color,legend_name] = get_color_scheme(color_scheme,files(cur_index));
                plot(files(cur_index).x+x_corr_factor,files(cur_index).y_cor_norm,lspec,'color',color,'linewidth',1,'DisplayName',char(legend_name));
                legend('-DynamicLegend');
                title(strcat(mutant_id,'C'),'FontName','arial','fontweight','normal')
            end
            return
    end
end

end



% function for determining the X offset required for alignment. 
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

% add styles here, and call them from the plot_mutant function. 
function plot_aesthetics(varargin)
switch varargin{1}
    case 'figure'
        %aesthetics:
        xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
        set(gca,'yticklabel','');
        axis([343 353 -inf inf])
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 1024, 800], 'PaperUnits', 'centimeters','PaperSize', [21, 29.7])
        set (gcf,'InvertHardcopy','off','Color',[1 1 1])
        grid off
        box off
    case 'default'
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 1024, 800], 'PaperUnits', 'centimeters','PaperSize', [21, 29.7])
        grid off
        set (gcf,'InvertHardcopy','off','Color',[1 1 1])
        axis([342 355 -inf inf])
        xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
    otherwise
       return
end
end


% add color schemes here, and call them within the plot_mutant function.
function [lspec,color,legend_name] = get_color_scheme(color_scheme,S)
identifier = strcat(S.state,S.antibody);
 
switch iscell(color_scheme)
    case 1
        lspec = '-';
        color = color_scheme{find(ismember(color_scheme,identifier))+1};
        legend_name = strcat(S.mutant,'C-',S.state,'-',S.antibody);
        if isempty(color)
            color = get_color_scheme('default',S);
        end
    
    otherwise
    switch color_scheme
        %add more cases as required
        case 'dualbinding'
            styles = {'monomerApo','-','monomer4b12','-','polymerApo','-','polymer5e3','-','monomer5e3','-','monomer4b12-5e3','-'};
            colors = {'monomerApo','k','monomer4b12','y','polymerApo','k','polymer5e3','r','monomer5e3','b','monomer4b12-5e3','r'};
            legend_name = strcat(S.mutant,'C-',S.state,'-',S.antibody);
        case '4b12figure'
            styles = {'monomerApo','-','monomer4b12','-','polymerApo','-','polymer5e3','-','monomer5e3','-','monomer4b12-5e3','-'};
            colors = {'monomerApo','k','monomer4b12','r','polymerApo','k','polymer5e3','r','monomer5e3','g','monomer4b12-5e3','b'};
            legend_name = S.antibody;
        case '5e3figure'
            styles = {'monomerApo','-','monomer4b12','-','polymerApo','-','polymer5e3','-','monomer5e3','-','monomer4b12-5e3','-'};
            colors = {'monomerApo','b','monomer4b12','r','polymerApo','k','polymer5e3','r','monomer5e3','','monomer4b12-5e3','-'};
            if strcmp(S.antibody,'Apo')
                legend_name = S.state;
            else
                legend_name = S.antibody;
            end
        case 'default'
            styles = {'monomerApo','-','monomer4b12','-','polymerApo','-','polymer5e3','-','monomer5e3','-','monomer4b12-5e3','-'};
            colors = {'monomerApo','k','monomer4b12','r','polymerApo','k','polymer5e3','r','monomer5e3','g','monomer4b12-5e3','b'};
            legend_name = strcat(S.mutant,'C-',S.state,'-',S.antibody);
        otherwise
            styles = {'monomerApo','-','monomer4b12','-','polymerApo','-','polymer5e3','-','monomer5e3','-','monomer4b12-5e3','-'};
            colors = {'monomerApo','k','monomer4b12','r','polymerApo','k','polymer5e3','r','monomer5e3','g','monomer4b12-5e3','b'};
            legend_name = strcat(S.mutant,'C-',S.state,'-',S.antibody);
    end
    lspec = styles{find(ismember(styles,identifier))+1};
    color = colors{find(ismember(colors,identifier))+1};
end
end

