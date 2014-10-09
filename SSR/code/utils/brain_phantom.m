function [tlobe,lobe,wm_mask,mask,tform] = brain_phantom( varargin )
%
% Return parameters:
%   tlobe     : deformed image
%   lobe      : original image w/o deformation
%   wm_mask   : white matter mask
%   mask      : brain mask
%   tform     : transformations for each tlobe image
%
% 'lobe' type images should be referenced in publications as: 
%    Internet brain segmentation repository (IBSR) [Online]. Available: http://www.cma.mgh.harvard.edu/ibsr.


Nt = 11;
gradient = 1;
GM = 80;
m = [256 256];
omega = [0 64 0 64]; 
rand_disp_optn={'p',[4 5]}; % spline parameters
lobe_type = 8;
parametric = 0;

% Parse varargin
for k=1:2:length( varargin ),
  eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end;

max_coeff=0.2/(128/m(1));

% Gradient
switch(gradient)
    case 0 % uniform
        WM=linspace(60,180,Nt-1)';
        WM=[WM(ceil(Nt/2)); WM];
    case 1 % increasing gradient
        WM(:,1)=linspace(200,0,Nt); % slopes
        WM(:,2)=ones(Nt,1).*50; % intercepts        
end

% Longitudinal displacement
rng('shuffle'); % random number generator is seeded with the same number when matlab starts, this reseeds it

if(parametric)
    wc_type = 'affine2D';
 
    max_disp = m(1)/20;
    aff_std  = 0.05;
    disp_std = 4;
    
    aff  = repmat([1 0 0 1]',1,Nt) + (randn(4,Nt)*aff_std);
    disp = randn(2,Nt)*disp_std*0;
    
    wc = [aff(1:2,:); disp(1,:); aff(3:4,:); disp(2,:)];
else
    wc_type = 'splineTransformation2D';

    wc=randn(prod(rand_disp_optn{2})*2,1);
    for ii=2:Nt
        wc(:,ii) = wc(:,ii-1) + randn(prod(rand_disp_optn{2})*2,1);
    end
    wc=wc./max(abs(wc(:))).*max_coeff;
end

[im,wm_mask,mask] = brainPhantom('type','lobe','lobe_type',lobe_type,'lobe_vent',1,'gradient',WM(1,:),'GMpoly',GM,'noiseImage',1,'snr',0.05,'m',m);        
lobe = im{1};
tlobe = im{1};

for pp=2:Nt
    % original image
    im = brainPhantom('type','lobe','lobe_type',lobe_type,'lobe_vent',1,'gradient',WM(pp,:),'GMpoly',GM,'noiseImage',1,'snr',0.05,'m',m);        
    lobe(:,:,pp) = im{1};

    % deform image 
    tform{pp} = struct( 'wc_type',wc_type,'wc',wc(:,pp),'yc_type',[],'yc',[],'omega',omega,'m',m );
    if( ~parametric )
        tform{pp}.params = rand_disp_optn;
    end
    
    im = brainPhantom('type','lobe','lobe_type',lobe_type,'lobe_vent',1,'gradient',WM(pp,:),'GMpoly',GM,'noiseImage',1,'snr',0.05,'m',m,'tform',tform{pp});        
    tlobe(:,:,pp) = im{1};
end

% reverse order
lobe  = lobe(:,:,end:-1:1);
tlobe = tlobe(:,:,end:-1:1);
tform = tform(end:-1:1);



