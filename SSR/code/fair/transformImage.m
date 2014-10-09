function [imt,yc,xc]=transformImage(im, tform, varargin)

% ---- Default Parameters ---- %

defaults.m              = size( im );
ndim                    = numel( defaults.m );
defaults.inter          = ['linearInter' num2str(ndim) 'D'];
defaults.regularizer    = 'moments';
defaults.theta          = 1e-1;
defaults.apply_wc       = 1;

if(isfield(tform,'omega'))
    defaults.omega = tform.omega;
else
    defaults.omega = [0 size(im,1) 0 size(im,2)];
end

% parse arguments
options = parseArgs( defaults, varargin );

if( isempty(tform) )
    imt=im;
    return;
end

% ---- Transformation ---- %

inter('reset','inter',options.inter,'regularizer',options.regularizer,'theta',options.theta);
trafo('reset','trafo',tform.wc_type,'omega',options.omega,'m',options.m);
if( isfield(tform,'params') )
    trafo('set',tform.params{:});
end

% resize image if necessary
im = imgresize( im, tform.omega, options.m );

% resize tform.yc if necessary
tform = tform_resize( tform, options.m );

xc = getCenteredGrid( tform.omega, options.m );
yc = xc;

% parametric transformation
if( isempty(tform.yc) ), options.apply_wc = 1; end
if( ~isempty(tform.wc) && options.apply_wc )
    yc = trafo(tform.wc, xc);   % compute the transformed points
end

% non-parametric deformation
if( isfield(tform,'yc') && ~isempty(tform.yc) )
    yc = yc + (center(tform.yc,tform.m) - xc); % 
end

imt = reshape(inter(im, options.omega, yc), options.m);

