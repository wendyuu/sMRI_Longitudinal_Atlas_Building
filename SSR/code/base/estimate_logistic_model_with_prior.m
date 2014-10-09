function [model,npoints] = estimate_logistic_model_with_prior( ti, yi, mask, varargin )
%function [model,npoints] = estimate_logistic_model_with_prior( ti, yi, mask, w, tip, yip, gammas, ymin, ymax, varargin )

%  Parse varargin
for k=1:2:length( varargin ) 
    eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end

m = size( mask );
model = zeros( [prod(m) 2] );
nt = numel( ti );

% stack images as 1D arrays
yi = reshape( yi, [], nt );
mask = reshape( mask, [], 1 );

if( exist('prior', 'var') )
    w = prior.weight;
    tip = prior.t;
    yip = reshape( prior.est, [], numel(prior.t) );
    gammas = reshape( prior.gammas, [], numel(prior.t) );
    ymin = prior.ymin;
    ymax = prior.ymax;
else
    w = 0;
    tip = 0;
    yip = zeros( size(yi,1), 1 );
    gammas = zeros( size(yi,1), 1 );
end

% normalize images to [0 1] -> c = 1
yi = (yi - ymin) / (ymax - ymin);
c = 1;

mask_idx = find( mask );
npoints = numel( mask_idx );

ct = 1;
tstart = tic;
for ii=mask_idx'
    [bhat, khat] = logistic_fit_with_prior( ti, yi(ii,:), c, w, tip, yip(ii,:), gammas(ii,:) );
    model(ii,:) = [bhat khat];
    if( ct == 100 )
        tend = toc(tstart);
        fprintf( 'First 100 logistic fits took %.2f secs. Estimated total time: %.2f min (%i points).\n', tend, npoints/6000*tend, npoints ); 
    end
    
    ct = ct + 1;
end

model = squeeze( reshape( model, [m 2] ) );