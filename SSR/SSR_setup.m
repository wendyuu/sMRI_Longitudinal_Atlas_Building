function SSR_setup

% Setup FAIR
FAIR_Dir = '~/research/apps/FAIR-Sandbox/FAIR';
run( fullfile(FAIR_Dir, 'startup.m') );

% Setup SSR
dirs = {'code','examples'};
for i=1:numel( dirs )
  if (exist(dirs{i}) == 7)
    path(path,genpath(dirs{i}));
  end
end
clear dirs;
