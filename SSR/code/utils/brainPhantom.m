function [im, wm_mask, mask, bmask, gm_mask] = brainPhantom( varargin )
rootdir='/home/icsapo/research/projects/LongitudinalAtlasBuilding/matlab/registration_test';

addpath(fullfile(rootdir,'utils'));

% ---- Default Parameters ---- %
wm_mask=[];
mask=[];
bmask=[];

% type
types = {'phantom','brainweb','lobe'};
defaults.type               = types{2};
% geometry
defaults.m                  = [100 100];                % image size
defaults.WMrratio           = [0.4 0.33];               % WM to brain radii ratios
% appearance
defaults.WMpoly             = 0;                        % WM intensity polynomial
defaults.GMpoly             = 0;                        % GM intensity polynomial
defaults.BG                 = 0;                        % background intensity
defaults.t                  = 0;                        % time
defaults.gradient           = [];

% noise
defaults.noiseImage         = 0;                        % add global Rician noise
defaults.noiseImageSigma    = 5;                        % sigma of global Rician noise
defaults.snr                = [];
defaults.noiseGeometric     = 0;                        % add geometric noise
defaults.noiseTissue        = 0;                        % add tissue intensity noise
defaults.noiseTissueSigma   = 5;                        % sigma of tissue intensity noise
% transformation
defaults.tform              = {};                       % transformation
% multilevel
defaults.multilevel         = 0;
defaults.minLevel           = 3;
% lobe phantom
defaults.lobe_type          = 4;
defaults.lobe_vent          = 1;
defaults.crop               = 16;
% brainweb phantom
defaults.brainweb_skull     = 0;
defaults.brainweb_pad       = 0;
defaults.brainweb_slices    = 89;
defaults.dimension          = 2;
defaults.norm               = 1;
defaults.normall            = 1;

% parse arguments
options = parseArgs( defaults, varargin );

if( options.type==3 & ~isfield(options, 'radii') )
    options.radii = ceil(options.m./[3 4]);   % brain radii
end

% if omega is not supplied, set to image size
if( ~isfield(options, 'omega') )
    switch( options.dimension)
        case 2, options.omega = [0 options.m(1) 0 options.m(2)];
        case 3, options.omega = [0 options.m(1) 0 options.m(2) 0 options.m(3)];
    end
end





% ---- Draw Brain Phantom ---- %

% GM and WM intensities
GM = polyval( options.GMpoly, options.t );
WM = polyval( options.WMpoly, options.t );

% Myelination simulation - we use similar contrast to monkey data
WMmin = 43; % from monkey data
WMmax = 116; % from monkey data
WMdiff = WMmax-WMmin;
GMnormratio = 1.7; % we need to divide the whole image by this to get the correct GM intensities (65)
if( isfield(options,'myelination') )
    WM = repmat(WMmin*1.7,1,numel(options.t));
end
% add tissue intensity noise                     
if(options.noiseTissue)
    GM = GM + ricianNoised(GM, options.noiseTissueSigma);
    WM = WM + ricianNoised(WM, options.noiseTissueSigma);
end

if( options.type==3)
    center = options.m./2;
    WMoffset = [0 options.radii(2)*0.35];
    WMradii  = ceil(options.radii.*options.WMrratio);
end

switch( options.type )
    case 'brainweb'
        t1o = imrotate(double(getfield(load_nii(fullfile(rootdir,'data/t1_icbm_normal_1mm_pn0_rf0.nii')),'img')),180);
        seg = imrotate(getfield(load_nii(fullfile(rootdir,'data/phantom_1.0mm_normal_crisp.nii')),'img'),180);

        if( options.dimension == 2 )
            t1o=t1o(:,:,options.brainweb_slices);
            seg=seg(:,:,options.brainweb_slices);
        
            if(options.brainweb_pad)
                pad=options.brainweb_pad;
                s=size(t1o);
                pt1o=zeros(s+2*pad);
                pseg=zeros(s+2*pad);
                pt1o(pad+1:s(1)+pad,pad+1:s(2)+pad)=t1o;
                pseg(pad+1:s(1)+pad,pad+1:s(2)+pad)=seg;
                t1o=pt1o;
                seg=pseg;
            end
        end
        
        % normalize to [0 255]
        t1o_min = min(t1o(:));
        t1o_max = double(max(t1o(:)));

        t1o = (t1o-t1o_min)./(t1o_max/255);

        if( ~options.brainweb_skull )
            mask = (seg==1 |  seg==2 | seg==3 | seg==8 );
            bmask = ( seg==2 | seg==3 | seg==8 );
            t1o = t1o.*mask;
        else
            mask = seg > 0;
            bmask = (seg>0 & seg~=1);
        end

        if( options.BG ~= 0 )
            t1o(~mask) = options.BG;
        end
        
        % tissue masks
        gm_mask = seg==2; % mean = 117
        wm_mask = seg==3; % mean = 155
        
        % resize images if necessary
        if( ~isequal(size(t1o), options.m) )
            if( options.dimension == 3 )
                rs = @(im) shiftdim(imresize(shiftdim(imresize(im,[options.m(1) options.m(2)],'nearest'),1),[options.m(2) options.m(3)],'nearest'),2);
                                
                t1o     = rs( t1o );      
                gm_mask = rs( gm_mask );      
                wm_mask = rs( wm_mask );  
                mask    = rs( mask );      
                bmask   = rs( bmask );      
            else
                t1o=imresize(t1o,options.m,'nearest');
                gm_mask=imresize(gm_mask,options.m,'nearest');
                wm_mask=imresize(wm_mask,options.m,'nearest');
                mask=imresize(mask,options.m,'nearest');
                bmask=imresize(bmask,options.m,'nearest');
            end
        end
        
        % tissue means
        gm_mean = double(mean(t1o(gm_mask)));
        wm_mean = double(mean(t1o(wm_mask)));
    case 'lobe'
        f=fopen(sprintf('%i.SNR0.img',options.lobe_type));
        if(options.lobe_type ~= 8)
            rot = 90;
        else
            rot = 180;
        end
        lobe=imrotate(double(reshape(fread(f,'uint16','b'),256,256)),rot);
        mask=lobe>0;
        if(options.lobe_vent)
            circ=getImplicitEllipse( zeros(256,256), [128 129], [15 15], 'fill', 1 );
            lobe(circ==1)=0;
        end
        bmask=lobe>0;
        % crop
        if(options.crop)
            lobe=lobe(options.crop:end-options.crop,options.crop:end-options.crop);
            mask=mask(options.crop:end-options.crop,options.crop:end-options.crop);
            bmask=bmask(options.crop:end-options.crop,options.crop:end-options.crop);
        end
        lobe=imresize(lobe,options.m,'nearest');
        mask=imresize(mask,options.m,'nearest');
        bmask=imresize(bmask,options.m,'nearest');
        wm_mask=lobe==254;
        gm_mask=lobe==128;
end

for ii=1:numel(options.t)
    switch( options.type )
        case 'phantom'
            % GM ellipse
            im{ii} = getImplicitEllipse( ones(options.m)*options.BG, center, options.radii, 'fill', GM(ii) );

            % add WM ellipses
            im{ii} = getImplicitEllipse( im{ii}, center+WMoffset, WMradii, 'fill', WM(ii), 'gradient', options.gradient );                      
            im{ii} = getImplicitEllipse( im{ii}, center-WMoffset, WMradii, 'fill', WM(ii), 'gradient', options.gradient );
        case 'brainweb'
           t1=t1o;
           blur_width=3;
           
           % gradient
            if( ~isempty(options.gradient) )
                   if( options.dimension == 2 )
                        bbox = regionprops(wm_mask,'BoundingBox');
                        if(isstruct(options.gradient)) % spline
                            sp = ppval(options.gradient,linspace(0,1,bbox(1).BoundingBox(4)));
                            if(options.norm)
                                grad = sp(end:-1:1)-mean(sp);
                                grad = grad./max(grad);
                                grad = grad.*options.norm(2) + mean(sp)/options.norm(1);
                            else
                                grad = sp(end:-1:1)/mean(sp);
                            end
                        else % polynomial
                            grad = polyval(options.gradient,linspace(0,1,bbox(1).BoundingBox(4)))/wm_mean;
                        end
                        xt = ceil(bbox(1).BoundingBox(2));
                        xb = xt + bbox(1).BoundingBox(4)-1;
                        wm_grad=zeros(size(wm_mask));
                        wm_grad(xt:xb,:) = double(wm_mask(xt:xb,:)).*repmat(grad',1,size(t1,2));
                        t1(wm_mask) = t1(wm_mask).*wm_grad(wm_mask);
                        t1b=imfilter(t1,fspecial('gaussian',3,2));
                        m=imdilate(wm_mask,ones(blur_width))-imerode(wm_mask,ones(blur_width));
                        t1(m>0)=t1b(m>0);
                   end
            elseif(WM(ii)>0)
                t1(wm_mask) = t1(wm_mask)-wm_mean+WM(ii);
                t1b=imfilter(t1,fspecial('gaussian',3,2));
                m=imdilate(wm_mask,ones(blur_width))-imerode(wm_mask,ones(blur_width));
                t1(m>0)=t1b(m>0);
            end
            
            % GM
            if(GM(ii)>0)
                t1(gm_mask) = t1(gm_mask)-gm_mean+GM(ii);
                t1b=imfilter(t1,fspecial('gaussian',3,2));
                m=imdilate(gm_mask,ones(blur_width))-imerode(gm_mask,ones(blur_width));
                t1(m>0)=t1b(m>0);
            end
            
            % Myelination
            if( isfield(options,'myelination') )
                t1=t1./GMnormratio;
                
                r = max(options.m)/2 * options.t(ii); % radius of myelination, depends on time (t=0 - no; t=1 - full)
                myelin = getImplicitEllipse( zeros(options.m), options.m./2, repmat(r,1,options.dimension)); 
                myelin = double(~myelin).*WMdiff;
                
                sigma = 10;
                h = zeros(1,1,10);
                h(1,1,:) = fspecial('gaussian', [1 10], sigma);

                copt = 'same';
                imfilter3D = @(im,h) convn(convn(convn(im,h,copt),shiftdim(h,1),copt),shiftdim(h,2),copt);

                myelinb = WMmax - imfilter3D( myelin, h );
                
                t1(wm_mask) = myelinb(wm_mask);
            end
            
            im{ii}=t1.*options.normall;
        case 'lobe'
            lo =lobe;
            
            % gradient
            if( ~isempty(options.gradient) )
                bbox = regionprops(wm_mask,'BoundingBox');
                if(isstruct(options.gradient))
                    % spline
                    grad = ppval(options.gradient,linspace(0,1,bbox(1).BoundingBox(4)+1));
                else
                    % polynomial
                    grad = polyval(options.gradient,linspace(0,1,bbox(1).BoundingBox(4)+1));
                end
                xt = floor(bbox(1).BoundingBox(2));
                xb = xt + bbox(1).BoundingBox(4);
                wm_grad=zeros(size(wm_mask));
                wm_grad(xt:xb,:) = double(wm_mask(xt:xb,:)).*repmat(grad',1,size(lo,2));
                lo(wm_mask) = wm_grad(wm_mask);
            elseif(WM(ii)>0)
                lo(wm_mask) = WM(ii);
            end
            
            % GM
            if(GM(ii)>0)
                lo(gm_mask) = GM(ii);
            end
            
            im{ii}=lo;
   end
            
    % convert to right-handed coordinate system (from left-handed MATLAB)
    im{ii} = imrotate(im{ii},-90);
    wm_mask = imrotate(wm_mask,-90);
    mask = imrotate(mask,-90);
    bmask = imrotate(bmask,-90);
    % apply transform
    if( ~isempty(options.tform) )
        im{ii} = transformImage(im{ii}, options.tform);
        wm_mask= transformImage(wm_mask, options.tform)>0.5;
        mask= transformImage(mask, options.tform)>0.5;
        bmask= transformImage(bmask, options.tform)>0.5;
   end

    % add global Rician noise                     
    if(options.noiseImage)
        if(~isempty(options.snr))
           sig=max(im{ii}(:))*options.snr;
        else
           sig= options.noiseImageSigma;
        end
        im{ii} = round(ricianNoised( im{ii}, sig ));
    end
    
    if(options.multilevel)
        im{ii} = getMultilevel(im(ii), options.omega, options.m, 'minLevel', options.minLevel, 'fig', 0);
    end
end


