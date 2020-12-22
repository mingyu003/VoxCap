function [res_diff]=finite_difference_coefficients(dx,num_diff_pnts,src_cen,obs_cen,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts)
%the function is used for 1 order partial derivative
h=0.1*dx;
if num_diff_pnts==2 % -1, 0, 1
    
    obs_cen_term1=obs_cen+unit_normal_obs*h;
    obs_cen_term2=obs_cen-unit_normal_obs*h;
    term1=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term1,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term2=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term2,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    res_diff=(term1-term2)/(2*h);
    
elseif num_diff_pnts==4 % -2, -1, 0, 1, 2
    
    obs_cen_term1=obs_cen-unit_normal_obs*2*h;
    obs_cen_term2=obs_cen-unit_normal_obs*h;
    obs_cen_term3=obs_cen+unit_normal_obs*h;
    obs_cen_term4=obs_cen+unit_normal_obs*2*h;    
    term1=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term1,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term2=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term2,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term3=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term3,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term4=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term4,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    res_diff=(term1-8*term2+8*term3-term4)/(12*h);
    
elseif num_diff_pnts==6 % -3, -2, -1, 0, 1, 2, 3
    
    obs_cen_term1=obs_cen-unit_normal_obs*3*h;
    obs_cen_term2=obs_cen-unit_normal_obs*2*h;
    obs_cen_term3=obs_cen-unit_normal_obs*h;
    obs_cen_term4=obs_cen+unit_normal_obs*h;
    obs_cen_term5=obs_cen+unit_normal_obs*2*h;
    obs_cen_term6=obs_cen+unit_normal_obs*3*h;
    
    term1=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term1,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term2=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term2,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term3=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term3,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term4=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term4,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term5=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term5,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    term6=compute_1overR_quad_mingyu(dx,src_cen,obs_cen_term6,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts);
    
    res_diff = (-term1+9*term2-45*term3+45*term4-9*term5+term6)/(60*h);
    
else
    error('Error! number of sample points should be 2 4 or 6')
end
