function [efield_hybrid] = compute_1over_R_partial_differentiation_integral(dx,src_cen,obs_cen,src_pnl_normal,obs_pnl_normal,smpl_pnts,smpl_wghts,num_diff_pnts,epsr,epsc)

% This routine calls the analytical code or quadrature code due to the
% distance between source and observer panels. Quadrature code gives
% accurate results for the far interactions & inaccurate results for the near 
% interactions. On the other hand, the analytical code gives accurate results for 
% the near interactions and inaccurate results for the far interactions due
% to phase canceling.
if length(src_pnl_normal)==1
    if (src_pnl_normal==1); src_pnl_normal=[1,0,0];
    elseif (src_pnl_normal == 2); src_pnl_normal=[0,1,0];
    elseif (src_pnl_normal == 3); src_pnl_normal=[0,0,1];
    elseif (src_pnl_normal == -1); src_pnl_normal=[-1,0,0];
    elseif (src_pnl_normal == -2); src_pnl_normal=[0,-1,0];
    elseif (src_pnl_normal == -3); src_pnl_normal=[0,0,-1];
    end
end
if length(obs_pnl_normal)==1
    if (obs_pnl_normal==1); obs_pnl_normal=[1,0,0];
    elseif (obs_pnl_normal == 2); obs_pnl_normal=[0,1,0];
    elseif (obs_pnl_normal == 3); obs_pnl_normal=[0,0,1];
    elseif (obs_pnl_normal == -1); obs_pnl_normal=[-1,0,0];
    elseif (obs_pnl_normal == -2); obs_pnl_normal=[0,-1,0];
    elseif (obs_pnl_normal == -3); obs_pnl_normal=[0,0,-1];
    end
end

if (norm(src_cen-obs_cen) > 6.5*dx) % use quadrature formula
    [efield_hybrid]=finite_difference_coefficients(dx,num_diff_pnts,src_cen,obs_cen,src_pnl_normal,obs_pnl_normal,smpl_pnts,smpl_wghts);
%       [efield_hybrid]=compute_1overR_int_quad_Efield(dx,src_cen,src_pnl_normal_vect,obs_cen,obs_pnl_normal_vect,smpl_pnts,smpl_wghts,6);
else % use analytical formula
    [efield_hybrid] = analytical_derivation_1overR(dx,src_cen,obs_cen,src_pnl_normal,obs_pnl_normal,epsr,epsc);    
%     [efield_hybrid] = compute_1overR_int_analy_Efield(dx,src_cen,src_pnl_normal_vect,obs_cen,obs_pnl_normal_vect,epsr,epsc);
end

