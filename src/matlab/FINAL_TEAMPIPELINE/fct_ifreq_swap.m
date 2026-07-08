% Author(s): Marti Geoffrey
% Epsztein Lab 2016

% This function calculates indexes based on another frequency f2 given indexes
% if1 based on frequency f1 with lower approximation.

function if2 = fct_ifreq_swap(if1, f1, f2)



if2 = floor((f2 / f1)*(if1 - 1)) + 1;


% f1 = 5;
% if1 = [1 2 3 4 5 6 7 8 9 10];
% f2 = 3.5;
% time1 = 0:1/f1:15;
% time2 = 0:1/f2:15;
% 
% time1(if1)
% time2(if2)

