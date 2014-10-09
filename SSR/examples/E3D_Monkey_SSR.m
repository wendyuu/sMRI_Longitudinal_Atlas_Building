function E3D_Monkey_SSR

load('E3D_Monkey_SSR.mat');

% Inputs:
%  syntheticTimeSeries - longitudinal image series
%  mask                - white matter mask of the last time-point
%  t                   - time-points
%  m                   - image resolution
%  omega               - physical image dimension
[imgs_it, tform, ests] = SSR_registration( monkeyTimeSeries, 'mask', mask, 't', t, ...
                                           'm', m, 'm_est', [64 64 64], 'omega', omega, ...
                                           'ymin', ymin, 'ymax', ymax, ...
                                           'verbose', 3 );

