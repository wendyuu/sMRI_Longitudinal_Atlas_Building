function [imgs_its, tform, ests] = SSR_registration( imgs, varargin )
% Model based similarity measure registration
%
% OPTIONS:
%  'use_model',<0|1>    : 0=regular, 1=model based
%  'm',<vector>         : image size used for registration
%  'out',<string>       : file name for saving tform

% If we get a file name, load it
if( ischar(imgs) )
    load( imgs );
end

sz   = size( imgs );

% Defaults
ntp         = sz(end);             % number of time-points
tfin        = ntp;                 % last time point index
m_est       = sz(1:end-1);         % original image size
m           = m_est;               % image size used for registration (images will be resized if m ~= m0)
ndim        = numel( m );          % number of dimensions
target      = tfin;                % use the last time point as the target image

est_tps  = 1:ntp;

% FAIR parameters
dists = {{'SSD'},{'NCC'},{'MI','nT',16,'nR',16}};   % available distance measures
switch( ndim )
    case 2,alphas  = [1000,1e-4,1e-4];              % 2D alpha corresponding to distance measure         
    case 3,alphas  = [1000,1e-6,1e-7];              % 3D alpha corresponding to distance measure 
end
dist        = 1;                                    % which distance measure to use
plot        = @myFAIRplots;                         % which function to use to plot progress
minLevel    = 4;                                    % MLdata lowest resolution level
parametric  = 0;                                    % parametric or deformable registration
preregister = 0;                                    % do a parametric preregistration for deformable registration

% mSM parameters
maxiter     = 10;        % maximum number of iterations
dSSR_tol    = 1e-3;     % maximum difference between SSRs for convergence
use_model   = 1;        % use longitudnal model (0 = regular registration)           
verbose     = 1;        % output progress

% output parameters
save_last_iter_only = 1;

%  Parse varargin
for k=1:2:length( varargin ) 
    eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end

if( ~exist('alpha','var') )
    alpha = alphas(dist);
end

n = prod( m );  % number of pixels/voxels


% Print parameters

if( use_model )
    fprintf( '======== Model based registration ' )
    if( exist('prior','var') ), fprintf( 'with prior ' ); end    
    fprintf( '========\n');
else
    fprintf( '======== FAIR registration ========\n');
end

fprintf( '    dist          = %s\n', dists{dist}{:} );
fprintf( '    alpha         = %f\n', alpha );
fprintf( '    maxiter       = %i\n', maxiter );
fprintf( '    m             = %s\n', mat2str(m) );
if( exist('m_est','var') ), fprintf( '    m_est         = %s\n', mat2str(m_est) ); end
fprintf( '    omega         = %s\n', mat2str(omega) );
fprintf( '    t             = %s\n', mat2str(t) );
if( exist('prior','var') )
    fprintf( '    prior\n' );
    fprintf( '       weight      = %f\n', prior.weight ); 
    fprintf( '       time-points = %s\n', int2str(prior.t) ); 
    fprintf( '       ymin        = %i\n', prior.ymin ); 
    fprintf( '       ymax        = %i\n', prior.ymax ); 
end
fprintf( '\n\n' );


% Dimension dependent parameters

switch( ndim )
    case 2
        img = @(im,tt) im(:,:,tt);
        viewImage( 'reset','viewImage','viewImage2D','colormap','gray(256)' );
    case 3
        img = @(im,tt) im(:,:,:,tt);
        viewImage( 'reset','viewImage','imgmontage','colormap','gray(256)','direction','-zyx' );
end

% FAIR options

interpolator = ['linearInter' num2str(ndim) 'D'];
transform    = ['affine' num2str(ndim) 'D'];

inter(       'reset','inter', interpolator                                            );
trafo(       'reset','trafo', transform                                               );
distance(    'reset','distance', dists{dist}{:}                                       );
regularizer( 'reset','regularizer','mfElastic','alpha',alpha,'mu',1,'lambda',0        );

optn = { 'maxIter', 200, ...
         'plotIter', 0, ...
         'plotMLiter', 0, ...
         'out', 0, ...
         'NPIR', @lBFGS, ...
         'NPIRobj', @NPIRobjFctn, ...
         'Plots', plot, ...
         'parametric', preregister };

% Resize images if necessary
if( ~isequal(m_est, m) )
    sz = size( imgs );
    tmp = zeros( [m_est sz(end)] );
    
    for tt=1:sz(end)
        switch( ndim )
            case 2, tmp(:,:,tt)   = imgresize( img(imgs,tt), omega, m_est );
            case 3, tmp(:,:,:,tt) = imgresize( img(imgs,tt), omega, m_est );
        end    
    end
    
    imgs = tmp;
   
    if( exist('mask','var') )
        mask = imgresize( mask, omega, m_est, 'method', 'nearest');
        
        % Update varargin
        varargin = replace_cell_value_pair( varargin, 'mask', mask );
    end
    
    if( exist('prior','var')  )
        sz = size( prior.est );
        tmp = zeros( [m_est sz(end)] );
        for tt=1:sz(end)
            switch( ndim )
                case 2, tmp(:,:,tt)   = imgresize( img(prior.est,tt), omega, m_est );
                case 3, tmp(:,:,:,tt) = imgresize( img(prior.est,tt), omega, m_est );
            end    
        end        
        prior.est = tmp;
        
       sz = size( prior.est );
       tmp = zeros( [m_est sz(end)] );
       for tt=1:sz(end)
            switch( ndim )
                case 2, tmp(:,:,tt)   = imgresize( img(prior.gammas,tt), omega, m_est );
                case 3, tmp(:,:,:,tt) = imgresize( img(prior.gammas,tt), omega, m_est );
            end    
        end        
        prior.gammas = tmp;
        
       % Update varargin
       varargin = replace_cell_value_pair( varargin, 'prior', prior );
    end
    
    clear tmp;
end


% Iterative Estimation and Registration

if( verbose>0 )
    fprintf( 'ITER:    ' );
    for i=1:tfin-1
        fprintf( '      du%i', i );
    end
    fprintf( ' |       du\n' );
end

imgs_it = imgs;

if( exist('prior', 'var') )
    args = {'prior', prior};
else
    args = {'ymin', ymin, 'ymax', ymax};
end

for iter=1:maxiter
    
    if( use_model )
        
        % Estimate model    
        if( verbose>1 ), fprintf( '%2i: estimating intensity model...\n', iter ); tstart = tic; end
        %[model,nPoints] = estimate_intensity_model( imgs_it, t, omega, varargin{:} );
        
        
        [model,nPoints] = estimate_logistic_model_with_prior( t, img(imgs_it,est_tps), mask, args{:} );
        
        % Save results
        if( exist('prefix','var') )
            save( [prefix '_model' num2str(iter) '.mat'], 'model', 'alpha', 'dist', 't', 'est_tps' );
        end
        
        if( verbose>1 ), fprintf( '  took %.2f mins (%i points)\n', toc(tstart)/60, nPoints ); end

    end
    
    
    % Register image to model

    for tt=1:tfin-1
        
        if( verbose>1 ), fprintf( '    registering time-point %i to %i... ', tt, tfin ); tstart = tic; end
               
        % Generate target images

        if( use_model )
            est = apply_logistic_model( img(imgs,target), model, mask, t(tt), args{:} );
        else
            est = img(imgs,target);
        end
        
        switch( ndim )
            case 2, 
                ests(:,:,tt,iter)   = est;               
            case 3, 
                ests(:,:,:,tt,iter) = est;
        end
        
        MLdata = getMultilevel( {img(imgs,tt),est}, omega, m,'minLevel', minLevel, 'fig', 0 );
                
        
        % Register

        if( parametric )
            [wc,his] = MLPIR( MLdata, optn{:},'verbose', (verbose>2) );
            yc = [];
        else
            [yc,wc,his] = MLIR( MLdata, optn{:}, 'verbose', (verbose>2) );
            yc = center( yc, m );
        end
        
        tform(tt,iter) = struct( 'wc_type', trafo('get','trafo'), 'wc', wc, 'yc_type', regularizer('get','regularizer'), 'yc', yc, 'omega', MLdata{end}.omega, 'm', MLdata{end}.m, 'dist', dists(dist), 'alpha', alpha );

        % Apply transform
        
        switch( ndim )
            case 2, 
                imgs_it(:,:,tt)   = transformImage( img(imgs,tt), tform(tt,iter), 'apply_wc', 0, 'omega', omega, 'inter', interpolator);               
                ssr(tt) = sum( reshape(est - imgs_it(:,:,tt),[],1).^2 )/n;
                
                imgs_its(:,:,tt,iter) = imgs_it(:,:,tt);
            case 3, 
                imgs_it(:,:,:,tt) = transformImage( img(imgs,tt), tform(tt,iter), 'apply_wc', 0, 'omega', omega, 'inter', interpolator);                
                ssr(tt) = sum( reshape(est - imgs_it(:,:,:,tt),[],1).^2 )/n;
                
                imgs_its(:,:,:,tt,iter) = imgs_it(:,:,:,tt);
        end
        
        if( verbose>1 ), fprintf( 'took %.2f mins\n', toc(tstart)/60 ); end
        
    end
    
    SSR(iter) = sum( ssr );
        
    if( verbose>0 )
        fprintf( '------------------ iter %2i ------------------\n      DU: ', iter );
        
        fprintf( '     SSR: ' );
        for i=1:tfin-1
            fprintf( '%8.2f ', ssr(i) );
        end
        fprintf( '| %8.4f\n', SSR(iter) );
    end
        
        
    if( ~use_model ) % regular registration
        break;    
    end
    
    % Check convergence
    
    if( iter > 1 )
               
        % Compute SSD between deformation fields of current and previous iterations
        if( parametric )
            dy = sum( ([tform(:,iter-1).wc]-[tform(:,iter).wc]).^2, 1 );
        else
            dy = sum( ([tform(:,iter-1).yc]-[tform(:,iter).yc]).^2, 1 )./n;
        end
        du = sum( dy );
        
        if( verbose>0 )
            for i=1:tfin-1
                fprintf( '%8.2f ', dy(i) );
            end
            fprintf( '| %8.2f\n', du );
        end
        
        % STOPPING CRITERIUM: check if SSR descreased between iterations
        dSSR = SSR(iter-1) - SSR(iter);
        
        if( dSSR < 0 || dSSR < dSSR_tol )
            break;
        end
        
    else
        
        if( verbose>0 ), fprintf( 'iter %2i:\n', iter ); end
          
    end
end

% Save results
if( exist('prefix','var') )
    % save tforms only for the last iteration 
    if( save_last_iter_only )
        tform = tform(:,end);
        
        switch( ndim )
            case 2, imgs_its = imgs_its(:,:,:,end);
            case 3, imgs_its = imgs_its(:,:,:,:,end);
        end
    end
    
    dist = dists{dist};
    save( [prefix '_tform.mat'], 'tform', 'alpha', 'dist', 't' );
    save( [prefix '_reg.mat'], 't', 'imgs_its' );
    if(use_model), save( [prefix '_est.mat'], 't', 'ests' ); end
end
end

%% -------------------------------------------- 
function cellArray = replace_cell_value_pair( cellArray, name, value )

    idx = [];
    for i=1:numel(cellArray)
        if( strcmp( cellArray{i}, name) )
           idx(end+1:end+2) = [i i+1];
        end
    end
    cellArray(idx) = [];

    cellArray(end+1:end+2) = {name, value};
end