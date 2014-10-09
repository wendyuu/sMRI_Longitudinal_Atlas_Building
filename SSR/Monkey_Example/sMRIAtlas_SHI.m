function sMRIAtlas_SHI
proc_path = '/primate/wendyuu_home/Local_Build'
addpath(fullfile(proc_path,'SSR'),fullfile(proc_path,'SSR/code'),fullfile(proc_path,'SSR/code/utils')
root_dir = '/Animal/primate/SanchezEmory/BrainDevYerkes/Longitudinal_sMRI_Atlas';
subj= {'RBf13','RBo14','RCh13','RDf13','RDw13','REw13','RFo14','RFv13','RHm13','RIm13','RJe13','RKr14','RLc14','RMe14','RNc14','RNe14','RNn14','RNt14','ROn14','RPv13','RPy13','RRm13','RRv13','RSb14','RSt14','RTm13','RUb14','RUf14','RUl13','RUu14','RVf14','RWm14','RWs14','RYl13','RYm14','RZu14'};
age={'2weeks','3months','6months','12months','18months'};

s.imres = [250 300 200]; %image resolution
s.omega = [ 68.3500   82.0200   54.6800];
for ss=1:numel(subj)
    fprintf('---- Subject %s ----\n',subj{ss});
    subj_dir = subj{ss};
    maskf = fullfile(root_dir,'WM_BM',sprintf('%s/%s_%s_T1_WM.nrrd', subj_dir,subj{ss},age{5}));
    if (exist(maskf) ~= 0) %last time point exists
        [temp,meta] = nrrdRead( maskf ) ;
        s.mask = double(temp);
        stt = 1; %index for TimeSeries 
        for tt=1:5
            fprintf('time-point %i\n',tt);
            imgf = fullfile(root_dir,'GM_Histmatched',sprintf('%s/%s_%s_T1_bias_masked_gmHistMatch.nrrd', subj_dir,subj{ss},age{tt}));
            if (exist(imgf)~= 0)  
                [temp,meta] = nrrdRead( imgf ) ;
                img  = double(temp);
                s.EmoryMonkeyTimeSeries(:,:,:,stt) = img;
                stt = stt + 1;
            end
        end
    else
        fprintf('No last time point for Subject %s\n',subj{ss});
    end
    s.tm_pt = stt - 1; % how many time points
    structf = fullfile(root_dir,'Matlab_Struct',sprintf(['%s/' ...
    '%s_EmoryMonkeyTimeSeries_with_%f_TimePoints.mat'], subj_dir,subj{ss},s.tm_pt));
    outdir =  fullfile(root_dir,'Matlab_Struct',sprintf(['%s/' ...
    '%s_E2D_SSR_with_%f_TimePoints.mat'], subj_dir,subj{ss},s.tm_pt));
    save(structf, '-struct', 's');
    printf(['saving ' ...
    '%s_EmoryMonkeyTimeSeries_with_%f_TimePoints.mat'],subj_dir);
    outfile = fullfile(outdir, ['E2D_Synthetic_SSR_cluster_' num2str(alpha) '.mat'] );
    alpha = 750
 %   SSR_registration(s.EmoryMonkeyTimeSeries, ...
 'mask', s.mask, ...
                  't', s.tm_pt, ...
                  'm', s.imres, ...
                  'omega', s.omega, ...
                  'alpha', alpha, ...
                  'out', outfile, ...
                  'plot', [], ...
                  'verbose', 3 );
end


[imgs_it] = SSR_registration( EmoryMonkeyTimeSeries, 'mask', mask, 't', tm_pt, ...
                              'm', imres, 'omega', omega, ...
                              'verbose', 3 );

figure;
for ii=1:Nt
    % original images
    subplot( 2, Nt, ii );    
        imagesc( imrotate(EmoryMonkeyTimeSeries(:,:,:,ii), 90), [0 255] ); colormap gray; axis off;
        
    % registered images
    subplot( 2, Nt, ii+Nt );
        imagesc( imrotate(imgs_it(:,:,:,ii), 90), [0 255] ); colormap gray; axis off;
end
