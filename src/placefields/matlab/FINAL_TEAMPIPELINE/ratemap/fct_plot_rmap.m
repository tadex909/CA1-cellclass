% Author(s): Marti Geoffrey
% Epsztein Lab 2019




function varargout = fct_plot_rmap(varargin)

%%%%colormap PP%%%%
CM = [   0.968627450980392         0.984313725490196                         1;
         0.870588235294118          0.92156862745098         0.968627450980392;
         0.776470588235294         0.858823529411765         0.937254901960784;
         0.619607843137255         0.792156862745098         0.882352941176471;
         0.419607843137255         0.682352941176471          0.83921568627451;
         0.258823529411765         0.572549019607843         0.776470588235294;
         0.129411764705882         0.443137254901961         0.709803921568627;
        0.0313725490196078         0.317647058823529         0.611764705882353;
        0.0313725490196078         0.188235294117647         0.419607843137255];


if strcmp(varargin{1}, 'genpar')
    varargout{1} = gen_params;
    return
end

[rmap, pf, vprm, save_path] = dealin(varargin);

isnans = sum(isnan(rmap(1).fr_cx), 2) == length(rmap(1).prm.xbin) - 1 - 2*rmap(1).prm.xbin_rem;

if any(isnans)
    C = find(isnans, 1, 'first') - 1;
else
    C = rmap(1).prm.nb_cond;
end

xbin = rmap(1).prm.xbin;
xbin_rem = rmap(1).prm.xbin_rem;
xbin([1:xbin_rem (end -(xbin_rem - 1)):end]) = [];
B = length(xbin) - 1;

xlabval = floor(linspace(1, B, 5));
crap = xlabval;
crap(1)=0;
%xlabstr = arrayfun(@mat2str, xbin(xlabval), 'uni', 0);
xlabstr = arrayfun(@mat2str, crap+10, 'uni', 0);

ispc = cell2mat(arrayfun(@(X) sum(sum(X.ispf_cx, 2)) > 0, pf(:), 'uni', 0));

if vprm.onlypc
    F = sum(ispc);
    idcel = find(ispc);
else
    F = numel(rmap);
    idcel = 1:F;
end

h = figure;
fct_fullscreen(gcf)

if ischar(save_path)
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    for gg = 1:F
        plot_ratemap_cel(gg);
        varargout{1} = h;
        set(h, 'Visible', 'off');
        fct_save_figure(h, [save_path filesep 'Cell_' num2str(gg) '_RateMap_IsNorm' num2str(vprm.norm)], 'eps');
        fct_save_figure(h, [save_path filesep 'Cell_' num2str(gg) '_RateMap_IsNorm' num2str(vprm.norm)], 'jpg');
        
        arrayfun(@delete,findall(0,'type','axes'))
        arrayfun(@delete,findall(0,'type','text'))
        arrayfun(@delete,findall(0,'tag','Colorbar'))
    end
    return
end

ind_cur = 1;
plot_ratemap_cel(ind_cur);
 varargout{1} = h;

uicontrol('Style', 'pushbutton', 'String', 'Previous',...
    'Units', 'Normalized', ...
    'Position', [0.45 0.95 0.05 0.05],...
    'Callback', @but_previous_c);

uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Units', 'Normalized', ...
    'Position', [0.50 0.95 0.05 0.05],...
    'Callback', @but_next_c);

    function but_previous_c(~, ~)
        cla
        if ind_cur  == 1
            ind_cur = F;
        else
            ind_cur = ind_cur - 1;
        end
        
        arrayfun(@delete,findall(0,'type','axes'))
        arrayfun(@delete,findall(0,'type','text'))
        arrayfun(@delete,findall(0,'tag','Colorbar'))
        plot_ratemap_cel(ind_cur);
    end

    function but_next_c(~, ~)
        cla
        if ind_cur < F
            ind_cur = ind_cur + 1;
        else
            ind_cur = 1;
        end
        
        arrayfun(@delete,findall(0,'type','axes'))
        arrayfun(@delete,findall(0,'type','text'))
        arrayfun(@delete,findall(0,'tag','Colorbar'))
        plot_ratemap_cel(ind_cur);
    end

    function plot_ratemap_cel(ind_cur)
        g = idcel(ind_cur);
        
        for c = 1:C
            idx = rmap(g).prm.idcond_t == c;
            
            fr_tx = rmap(g).fr_s_tx(idx, :);
            fr_x = rmap(g).fr_s_cx(c, :);
            
            ispf_tx = pf(g).ispf_tx(idx, :);
            ispf_x = pf(g).ispf_cx(c, :);
            
            nbpf = max(ispf_x);
            if vprm.onlyonepf && (nbpf > 1)
                nbpf = 1;
            end
            nblap = sum(idx);
            
            base_line = mean(fr_x(fr_x <= prctile(fr_x, 10)));
            pf_line = 0.2*(max(fr_x) - base_line) + base_line;
            
            subplot(ceil(C / 2), 4, 2*(c-1) + 1)
            if vprm.norm
                fr_tx = fct_matnorm(fr_tx);
            end
            imagesc(fr_tx)
            
            % surf(gca, 1:B, 1:nblap, fr_tx)
            % shading interp;
            % axis tight;
            % view(0,90);
            
            hold on
            %colormap jet
            %colormap (CM)
            if c <3
            %colormap (vprm.colormapbnl)
            colormap(gca, vprm.colormapbnl);
            else
            %colormap (vprm.colormaponl)
            colormap(gca, vprm.colormaponl);
            end
            
            tt = colorbar;
            cpos = get(tt, 'Position');
            cpos(3) = 0.5*cpos(3);
            cpos(1) = cpos(1) + 0.03;
            set(tt, 'Position', cpos);
            
            %%%%PF
            for p = 1:nbpf
                idx = find(ispf_x == p);
                istart = idx(1);
                istop = idx(end);
                line([istart istart], [1 nblap],  'color', 'k', 'linewidth', 2, 'linestyle', '--')
                line([istop istop], [1 nblap],  'color', 'k', 'linewidth', 2, 'linestyle', '--')
                %%%individual PF
                %{
                for l = 1:nblap
                    idx = find(ispf_tx(l, :) == p);
                    if any(idx)
                        istart = idx(1);
                        istop = idx(end);
                        line([istart istart], [l-0.5 l+0.5],  'color', 'k', 'linewidth', 2, 'linestyle', '-')
                        line([istop istop], [l-0.5 l+0.5],  'color', 'k', 'linewidth', 2, 'linestyle', '-')
                    end
                end
                %}
                %%%%%%
            end
            
            xlabel('Spatial Bin')
            ylabel('Temporal Bin')
%             set(gca, 'xtick', xlabval, 'xticklabel', xlabstr) OLD GOOD FOR 0-80
            set(gca, 'xtick', xlabval, 'xticklabel', xlabval)
            info_txt = cell(1, 1);
            info_txt{1} = ['Cell ' num2str(g) '   Active:' num2str(rmap(g).feat(c).isactive) '   Fr: ' num2str(rmap(g).feat(c).fr, '%.1f')];
            info_txt{2} = ['Stb OddEven: ' num2str(rmap(g).feat(c).stb_oddeven, '%.2f') '   Mean: ' num2str(rmap(g).feat(c).stb_meancorr, '%.2f') '   AllPairs: ' num2str(rmap(g).feat(c).stb_allpairs, '%.2f')];
            info_txt{3} = ['SI: ' num2str(rmap(g).feat(c).si, '%.2f') '   MeanLap: ' num2str(rmap(g).feat(c).si_meanlap, '%.2f')];
            info_txt{4} = ['Sparsity: ' num2str(rmap(g).feat(c).sparsity, '%.2f')];
            
            title(info_txt)
                      
            subplot(ceil(C / 2), 4, 2*(c-1) + 2);
            
            plot(fr_x)
            hold on
                      
            plot(base_line*ones(1, B), 'r')
            plot(pf_line*ones(1, B), 'k')
            
%             Délimitations zones objets SPO
%             ylim('auto');
%             p = [0 5.75 15.25 21.75 35.5 42 54.25 60.75 74.25 80];%Coordonnées extrémités zone objet = 15cm au total
%             plot([p(1) p(1)], ylim , 'linestyle', '--', 'color', '[0.8500 0.3250 0.0980]');%Origami
%             plot([p(2) p(2)], ylim , 'linestyle', '--', 'color', '[0.8500 0.3250 0.0980]');
%             plot([p(3) p(3)], ylim , 'linestyle', '--', 'color', '[0.9290 0.6940 0.1250]');%Cube
%             plot([p(4) p(4)], ylim , 'linestyle', '--', 'color', '[0.9290 0.6940 0.1250]');
%             plot([p(5) p(5)], ylim , 'linestyle', '--', 'color', '[0.4940 0.1840 0.5560]');%Sablier
%             plot([p(6) p(6)], ylim , 'linestyle', '--', 'color', '[0.4940 0.1840 0.5560]');
%             plot([p(7) p(7)], ylim , 'linestyle', '--', 'color', '[0.4660 0.6740 0.1880]');%Maison
%             plot([p(8) p(8)], ylim , 'linestyle', '--', 'color', '[0.4660 0.6740 0.1880]');
%             plot([p(9) p(9)], ylim , 'linestyle', '--', 'color', '[0.3010 0.7450 0.9330]');%Sapin
%             plot([p(10) p(10)], ylim , 'linestyle', '--', 'color', '[0.3010 0.7450 0.9330]');
                       
            for p = 1:nbpf
                idx = find(ispf_x == p);
                if p == 1
                    plot(idx, fr_x(idx), 'or', 'markerfacecolor', 'r')
                else
                    plot(idx, fr_x(idx), 'ob', 'markerfacecolor', 'b')
                end
            end
            
            xlabel('Spatial Bin')
            ylabel('Firing Rate (Hz)')
            set(gca, 'xtick', xlabval, 'xticklabel', xlabstr)
            info_txt = cell(1, 1);
            if nbpf >= 1
                info_txt{1} = ['RateOut/In Field: ' num2str(pf(g).feat(c).outovin(1), '%.2f')];
                info_txt{2} = ['PF Width: ' num2str(pf(g).feat(c).size(1), '%.2f')];
                info_txt{3} = ['Mean PF/Lap Width: ' num2str(pf(g).feat(c).lap_size(1), '%.2f')];
                info_txt{4} = ['PF Deviation: ' num2str(pf(g).feat(c).lap_devwithmean(1), '%.2f')];
            else
                info_txt{1} = 'No Place Field';
            end
            
            title(info_txt)
           
        end
        
        set(gcf, 'color', 'w')
        
        fct_fullscreen(gcf)        
    end
end


function [rmap, pf, vprm, save_path] = dealin(X)

rmap = X{1};
pf = X{2};

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
loadcolmathSFN
load colormaprom.mat
vprm.colormapbnl = colormapbnl;
vprm.colormaponl = colormaponl;
vprm.norm = true;
vprm.onlypc = false;
vprm.onlyonepf = false;
end
