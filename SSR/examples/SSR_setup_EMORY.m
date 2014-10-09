function SSR_setup_EMORY

% Setup FAIR
FAIR_Dir = '/lustre/scr/w/e/wendyuu/FAIR';
run( fullfile(FAIR_Dir, 'startup.m') );

% Setup SSR
SSR_Dir = '/lustre/scr/w/e/wendyuu/SSR';
dirs = {'code','examples'};
for i=1:numel( dirs )
  if (exist(fullfile(SSR_Dir,dirs{i})) == 7)
    path(path,genpath(fullfile(SSR_Dir,dirs{i})))
  end
end
clear dirs;
