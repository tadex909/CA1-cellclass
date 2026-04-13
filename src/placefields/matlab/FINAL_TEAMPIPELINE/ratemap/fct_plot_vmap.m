% -----------------------------
% Written by MARTI Geoffrey
% 6/15
% 12/15
% 08/16
% -----------------------------



function varargout = fct_plot_vmap(varargin)

if strcmp(varargin{1}, 'genpar')
    varargout{1} = gen_params;
    return
end

[vmap, vprm, save_path] = dealin(varargin);

switch vprm.med_or_std
    case 1
        vcx = vmap.vmed_s_cx;
        vtx = vmap.vmed_s_tx;
    case 2
        vcx = vmap.vstd_s_cx;
        vtx = vmap.vstd_s_tx;
end

isnans = sum(isnan(vmap.vmed_cx), 2) == length(vmap.prm.xbin) - 1 - 2*vmap.prm.xbin_rem;

if any(isnans)
    C = find(isnans, 1, 'first') - 1;
else
    C = vmap.prm.nb_cond;
end

xbin = vmap.prm.xbin;
xbin_rem = vmap.prm.xbin_rem;
xbin([1:xbin_rem (end -(xbin_rem - 1)):end]) = [];
B = length(xbin) - 1;

xlabval = floor(linspace(1, B, 5));
xlabstr = arrayfun(@mat2str, xbin(xlabval), 'uni', 0);



h = figure;
fct_fullscreen(h)

if ischar(save_path)
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    
    plot_ratemap_cel;
    varargout{1} = h;
    set(h, 'Visible', 'off');
    fct_save_figure(h, [save_path filesep '_Vmap_IsNorm' num2str(vprm.norm)], 'jpg');
    
    arrayfun(@delete,findall(0,'type','axes'))
    arrayfun(@delete,findall(0,'type','text'))
    arrayfun(@delete,findall(0,'tag','Colorbar'))
    return
end


plot_ratemap_cel;
varargout{1} = h;

    function plot_ratemap_cel
        for c = 1:C
            idx = vmap.prm.idcond_t == c;
            
            fr_tx = vtx(idx, :);
            fr_x = vcx(c, :);
            
            
            %             nblap = sum(idx);
            
            base_line = mean(fr_x(fr_x <= prctile(fr_x, 10)));
            pf_line = 0.2*(max(fr_x) - base_line) + base_line;
            
            subplot(ceil(C / 2), 4, 2*(c-1) + 1)
            if vprm.norm
                fr_tx = fct_matnorm(fr_tx, 2, 'minmaxnorm');
            else
                fr_tx = (fr_tx - nanmin(fr_tx(:))) / (nanmax(fr_tx(:)) - nanmin(fr_tx(:)));
            end
            imagesc(fr_tx)
            
            hold on
            colormap jet
            
            tt = colorbar;
            cpos = get(tt, 'Position');
            cpos(3) = 0.5*cpos(3);
            cpos(1) = cpos(1) + 0.03;
            set(tt, 'Position', cpos);
            
            xlabel('Spatial Bin')
            ylabel('Temporal Bin')
            set(gca, 'xtick', xlabval, 'xticklabel', xlabstr)
            
            subplot(ceil(C / 2), 4, 2*(c-1) + 2)
            
            plot(fr_x)
            hold on
            
            
            plot(base_line*ones(1, B), 'r')
            plot(pf_line*ones(1, B), 'k')
            
        
            
            
            xlabel('Spatial Bin')
            ylabel(vprm.ylab)
            set(gca, 'xtick', xlabval, 'xticklabel', xlabstr)
        end
        
        
        set(gcf, 'color', 'w')
        
        fct_fullscreen(gcf)
        
        
    end

end


function [vmap, vprm, save_path] = dealin(X)

vmap = X{1};

if any(strcmp(X, 'vprm'))
    vprm = X{find(strcmp(X, 'vprm')) + 1};
else
    vprm = gen_params;
end

if any(strcmp(X, 'save'))
    save_path = X{find(strcmp(X, 'save')) + 1};
else
    save_path = [];
end

end

function vprm = gen_params
vprm.norm = false;
vprm.ylab = 'MeanV';
vprm.med_or_std = 1;
end




