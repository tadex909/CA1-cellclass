% Author(s): Marti Geoffrey
% Epsztein Lab 2019

function [] = fct_save_figure(hh, file_path_full, varargin)


set(hh, 'PaperPositionMode', 'auto')


folder_path = fileparts(file_path_full);


if ~(exist(folder_path, 'dir') == 7)
    mkdir(folder_path);
    disp(['The save folder ''' folder_path '''  has been created.'])
end


if nargin > 2
    for k = 1:numel(varargin)
        switch varargin{k}
            case 'jpg'
                if any(strcmp(varargin, 'poor'))
                    print( hh, '-djpeg', file_path_full, '-r100')
                else
                    print(hh, '-djpeg', file_path_full, '-r450')
%                     print(hh, '-djpeg', file_path_full, '-r2000')
                end
            case 'eps'
                print(hh, '-depsc2', '-painters', file_path_full)
            case 'fig'
                savefig(hh, file_path_full)
        end
    end
end


