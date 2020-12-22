% use recent DIRECTFN (faster)
use_recent_DIRECTFN=0;

% find the current folder
currentfolder = pwd;

% get system-dependent path separator
psep = pathsep();

% include all the needed path (addpath() adds recursively),
% but do NOT include the Git directory


p1  = genpath([currentfolder, filesep, 'circulant_tensor']);
p2  = genpath([currentfolder, filesep, 'one_over_R']);
p3  = genpath([currentfolder, filesep, 'src_iterative_solver']);
p4  = genpath([currentfolder, filesep, 'src_post_process']);
p5  = genpath([currentfolder, filesep, 'src_pre_process']);
p6  = genpath([currentfolder, filesep, 'Cross_Tuck_Toolbox']);
p9  = genpath([currentfolder, filesep, 'SVD_Tuck']);
p8  = genpath([currentfolder, filesep, 'preconditioner']);
p10  = genpath([currentfolder, filesep, 'prestored_Toeplitz']);
addpath([p1, psep,p2, psep, p3, psep, p4, psep, p5, psep, p6, psep, p8, psep , p9, psep, p10, psep]);

if (ispc)
  p7 = genpath([currentfolder, filesep, 'windows']);
else
  p7 = genpath([currentfolder, filesep, 'linux']);
end
addpath(p7);

% get system-dependent file separator
filesep = filesep();

currentfolder = pwd;


  

    
  
