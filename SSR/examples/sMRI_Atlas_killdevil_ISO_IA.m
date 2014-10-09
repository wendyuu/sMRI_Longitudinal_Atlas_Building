function sMRI_Atlas_killdevil_ISO_IA
proc_path = '/lustre/scr/w/e/wendyuu';
addpath(fullfile(proc_path,'SSR'),fullfile(proc_path,'SSR/code'),fullfile(proc_path,'SSR/code/utils'));
root_dir = '/lustre/scr/w/e/wendyuu/Longitudinal_sMRI_Atlas';
subj= {'RBf13','RBo14','RCh13','RDf13','RDw13','REw13','RFo14','RFv13','RHm13','RIm13','RJe13','RKr14','RLc14','RMe14','RNc14','RNe14','RNn14','RNt14','ROn14','RPv13','RPy13','RRm13','RRv13','RSb14','RSt14','RTm13','RUb14','RUf14','RUl13','RUu14','RVf14','RWm14','RWs14','RYl13','RYm14','RZu14'};
age={'2weeks','3months','6months','12months','18months'};
ages = {14,90,180,360,540};

s.imres = [256 256 256]; %image resolution
s.omega = [0 68.3500 0 82.0200 0 54.6800];
f = fopen( fullfile(proc_path,'SSR/examples/sMRI_Atlas_bsub_ISO_IA.csh'), 'w' );
fprintf( f, '# This file was gererated by sMRI_Atlas_killdevil_ISO_IA.m\n' ...
         );

for ss=1:numel(subj)
    fprintf('---- Subject %s ----\n',subj{ss});
    subj_dir = subj{ss};
    for tmp_t = 5 : -1 : 1
        maskf = fullfile(root_dir,'WM_BM_ISO',sprintf('%s/%s_%s_T1_WM_256iso.nrrd', subj_dir,subj{ss},age{tmp_t}));
        last_tm = age{tmp_t};
        if (exist(maskf) ~= 0) %last time point exists
            break;
        end
    end
    [temp,meta] = nrrdRead( maskf ) ;
    s.mask = double(temp);
    s.tm_pt = [];
    stt = 1; %index for TimeSeries 
    for tt=1:5
        fprintf('time-point %i\n',tt);
        imgf = fullfile(root_dir,'GM_Histmatched_ISO_IntensityAdjusted',sprintf('%s/%s_%s_T1_bias_masked_gmHistMatch_256iso_intensity_adjusted.nrrd', subj_dir,subj{ss},age{tt}));
        if (exist(imgf)~= 0)  
            [temp,meta] = nrrdRead( imgf ) ;
            img  = double(temp);
            s.EmoryMonkeyTimeSeries(:,:,:,stt) = img;
            s.tm_pt = [s.tm_pt,ages{tt}];
            stt = stt + 1;
        end
    end
    % how many time points
    structf = fullfile(root_dir,'Matlab_Struct_ISO_IntensityAdjusted',sprintf('%s/%s_EmoryMonkeyTimeSeries_with_%s_TimePoints_LastTime_%s_256iso_IntensityAdjusted.mat', subj_dir,subj{ss},num2str(stt-1),num2str(last_tm)));
    if (exist(structf) == 0 )
      save(structf, '-struct', 's');
      sprintf(['saving ' ...
      '%s_EmoryMonkeyTimeSeries_with_%s_TimePoints_256iso_IntensityAdjusted.mat'],subj_dir,num2str(stt-1)) 
    end
    fprintf(f,'bsub -M 10 -q week matlab -nojvm -nodisplay -singleCompThread -r "Run_SSR_EmoryMonkey(\'%s\')" -logfile %s\n',structf,strrep(structf,'.mat','_output.log'));
end