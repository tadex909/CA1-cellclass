% Author(s): Marti Geoffrey
% Epsztein Lab 2019

function [] = fct_fullscreen(hh)


set(hh, 'Units', 'normalized','Position',[0 0 1 1])

% Another way
% s = get(0, 'ScreenSize');
% set(hh, 'Position', [0 0 s(3) s(4)]);