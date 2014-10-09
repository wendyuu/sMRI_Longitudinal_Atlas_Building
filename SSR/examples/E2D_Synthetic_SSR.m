function E2D_Synthetic_SSR

load('E2D_Synthetic_SSR.mat');

% Inputs:
%  syntheticTimeSeries - longitudinal image series
%  mask                - white matter mask of the last time-point
%  t                   - time-points
%  m                   - image resolution
%  omega               - physical image dimension
[imgs_it, tform, ests] = SSR_registration( syntheticTimeSeries, 'mask', mask, 't', t, ...
                                           'm', m, 'm_est', [128 128], 'omega', omega, ...
                                           'ymin', 56, 'ymax', 254, ...
                                           'verbose', 3 );

figure;
for ii=1:Nt
    % original images
    subplot( 3, Nt, ii );    
        imagesc( imrotate(syntheticTimeSeries(:,:,ii), 90), [0 255] ); colormap gray; axis off;
end

for ii=1:Nt-1
    % estimated images
    subplot( 3, Nt, ii+Nt );    
        imagesc( imrotate(ests(:,:,ii,end), 90), [0 255] ); colormap gray; axis off;

    % registered images
    ycc = reshape( tform(end).yc, [tform(end).m 2] );
    subplot( 3, Nt, ii+2*Nt );
        imagesc( imrotate(imgs_it(:,:,ii,end), 90), [0 255] ); colormap gray; axis off;
        hold on;
            contour( ycc(:,:,1), 20, 'r', 'LineWidth', 1 );
            contour( ycc(:,:,2), 20, 'r', 'LineWidth', 1 );
        hold off;
end
