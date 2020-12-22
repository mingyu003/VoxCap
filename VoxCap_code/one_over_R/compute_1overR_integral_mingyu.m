function int_res=compute_1overR_integral_mingyu(dx,src_pnl_loc,obs_pnl_loc,src_pnl_normal,obs_pnl_normal,smpl_pnts,smpl_wghts)
%---------------------------------------------------------
%generate 1/R integral matrix
%--------------------------------------------------------
%%INPUTS
%dx: resolution namely the dimension of every voxel
%num_smpl: the number of integral sample points
%src_pnl_loc: center coordinates of source panel
%obs_pnl_loc: center coordinates of observer panel
%src_pnl_normal: unit normal vector of source panel
%obs_pnl_normal: unit normal vector of observer panel
%smpl_pnts,smpl_wghts: gauss points and weights
%%OUTPUT
%int_res: the result of i/R integral
%--------------------------------------------------------
if length(src_pnl_normal)==1
    if (abs(src_pnl_normal)==1); src_pnl_normal=[1,0,0];
    elseif (abs(src_pnl_normal) == 2); src_pnl_normal=[0,1,0];
    elseif (abs(src_pnl_normal) == 3); src_pnl_normal=[0,0,1];
    end
end
if length(obs_pnl_normal)==1
    if (abs(obs_pnl_normal)==1); obs_pnl_normal=[1,0,0];
    elseif (abs(obs_pnl_normal) == 2); obs_pnl_normal=[0,1,0];
    elseif (abs(obs_pnl_normal) == 3); obs_pnl_normal=[0,0,1];
    end
end
% analytic formula cannot be avaiable when the distance of two panels is
% large and gauss-legendre quadrature doesn't woek when  the distance of
% two panels is very small.
if (norm(obs_pnl_loc-src_pnl_loc) > 100*dx) % use quadrature formula
    int_res=compute_1overR_quad_mingyu(dx,src_pnl_loc,obs_pnl_loc,src_pnl_normal,obs_pnl_normal, smpl_pnts,smpl_wghts);
else % use analytical formula
    int_res=compute_1overR_analy_mingyu(dx,src_pnl_loc,obs_pnl_loc,src_pnl_normal,obs_pnl_normal);
end