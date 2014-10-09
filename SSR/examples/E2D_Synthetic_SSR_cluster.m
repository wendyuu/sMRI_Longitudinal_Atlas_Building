function E2D_Synthetic_SSR_cluster( index )
% Helper function for submitting 3 registration jobs (each with a different alpha)
% on the cluster (killdevil.unc.edu) using bsub.

outdir = '/nas02/home/i/c/icsapo';

load('E2D_Synthetic_SSR.mat');

alphas  = [ 500, 750, 1000 ];
alpha   = alphas(index);
outfile = fullfile( outdir, ['E2D_Synthetic_SSR_cluster_' num2str(alpha) '.mat'] );

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
SSR_registration( syntheticTimeSeries, ...
                  'mask', mask, ...
                  't', t, ...
                  'm', m, ...
                  'omega', omega, ...
                  'alpha', alpha, ...
                  'out', outfile, ...
                  'plot', [], ...
                  'verbose', 3 );

