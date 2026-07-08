function [theta, rbar, delta, sygma] = circmeanPP( DATA )
% [THETA, RBAR, DELTA, SYGMA] = CIRCMEAN( DATA ) 
%
% Function returns the non-weighted circular mean (THETA) of all of
% the elements in the vector DATA. Optionally returns the mean
% resultant length (rbar) and circular dispersion (delta) and circular standard error.
% from Fisher's book, equations origins is displayed in code
% PPLS 
%
if nargin == 0
  help circmean;
  return
end

DATA = DATA(:);
N = length(DATA);

C = sum( cos(DATA) );
S = sum( sin(DATA) );

theta = atan(S/C);
if C < 0
  theta = theta + pi;
elseif S < 0
  theta = atan(S/C) + 2*pi;% fisher p 31 eq. 2.9
end

if nargout >= 2
  rbar = sqrt( C^2 + S^2 )/N; % eq 2.7 and 2.10
end
if nargout >= 3
  moment2 = 1/N * sum( cos( 2 * (DATA - theta) ) ); %fisher's book p34 eq 2.27
  delta = (1 - moment2)/(2 * rbar^2);%fisher's book p34 eq 2.28
end

if nargout >=4
    sygma=sqrt(delta/N);
end;
