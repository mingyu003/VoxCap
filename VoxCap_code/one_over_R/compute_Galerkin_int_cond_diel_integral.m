function [Efield_or_pot] = compute_Galerkin_int_cond_diel_integral(dx,cen_src,norm_src,cen_obs,norm_obs,diel_or_cond,eps_inner,eps_outter,smpl_pnts,smpl_wghts,num_diff_pnts)
% diel_or_cond=0, the panel is corresponding conductor and we will
% calculate potential.
% diel_or_cond=1, the panel is corresponding dielectric and we will
% calculate electrical field with differentiation.
len=length(smpl_wghts);
if(diel_or_cond==0) % observer panel is conductor
    
    [Efield_or_pot] = compute_1overR_integral_mingyu(dx,cen_src,cen_obs,norm_src,norm_obs,smpl_pnts,smpl_wghts);
%     [Efield_or_pot] = compute_1overR_integral_mexn(dx,cen_src,cen_obs,norm_src,norm_obs,smpl_pnts,smpl_wghts,len);
    
elseif(diel_or_cond>=1) % observer panel is dielectric

    [Efield_or_pot] = compute_1over_R_partial_differentiation_integral(dx,cen_src,cen_obs,norm_src,norm_obs,smpl_pnts,smpl_wghts,num_diff_pnts,eps_inner,eps_outter);
%     [Efield_or_pot] =compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,norm_src,norm_obs,smpl_pnts,smpl_wghts,...
%                         len,num_diff_pnts);
else
    error('Invalid panel type - panel should be a conductor or dielectric')
    
end


end