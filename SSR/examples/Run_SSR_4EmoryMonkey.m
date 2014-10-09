function Run_SSR_4EmoryMonkey(f_EmoryMonkeyTimeSeries)

load(f_EmoryMonkeyTimeSeries);

% Inputs:
%  syntheticTimeSeries - longitudinal image series
%  mask                - white matter mask of the last time-point
%  t                   - time-points
%  m                   - image resolution
%  omega               - physical image dimension
%  alpha               - regularization weigth
%  outfile             - path to save results
%  plot                - [] = turn off plotting
%  verbose             - verbose level


alpha = 750;
outfile = strrep(f_EmoryMonkeyTimeSeries, '.mat', '_SSR_Output.mat');

SSR_registration( EmoryMonkeyTimeSeries, ...
                  'mask', mask, ...
                  't', tm_pt, ...
                  'm', imres, ...
                  'omega', omega, ...
                  'alpha', alpha, ...
                  'out', outfile, ...
                  'plot', [], ...
                  'verbose', 3 );

