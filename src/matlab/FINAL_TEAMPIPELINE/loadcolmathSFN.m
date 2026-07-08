% Optional colormap file. Use it when available, but do not fail if missing.
this_dir = fileparts(mfilename('fullpath'));
candidate_paths = {
    fullfile(this_dir, 'colormaprom.mat')
    fullfile(pwd, 'colormaprom.mat')
    };

loaded_colormap = false;
for ii = 1:numel(candidate_paths)
    if exist(candidate_paths{ii}, 'file') == 2
        load(candidate_paths{ii});
        loaded_colormap = true;
        break
    end
end

if ~loaded_colormap
    warning('loadcolmathSFN:MissingColormapFile', ...
        'colormaprom.mat not found. Continuing with built-in color_set/color_patch only.');
end

color_set = [65 97 120 ; 120 175 218 ; 251 183 86 ; 127 55 169 ; 224 176 255] / 255;

color_patch = [[180 20 1] ; [123 125 132] ; [180 20 1] ; [123 125 132]]/ 255;