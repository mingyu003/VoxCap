function [G_block] = obtain_z_mat_block(dx,ind_arr_obs,ind_arr_src,geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts)

%inputs:
%ind_arr_obs,ind_arr_src: dimensions of z_mat
%dx: resolution
%smpl_pnts,smpl_wghts: Gauss points and weights
num_src=length(ind_arr_src);
num_obs=length(ind_arr_obs);
G_block(:,:)=zeros(num_obs,num_src);
parfor kk=1:num_obs % observer panel
    ind_obs=ind_arr_obs(kk);
    for ll=1:num_src% source panel
        ind_src=ind_arr_src(ll);
        [G_block(kk,ll)]=compute_Galerkin_int_cond_diel_integral(dx, ...
            geom_bndry_panels(ind_src,1:3),geom_bndry_panels(ind_src,4), ...
            geom_bndry_panels(ind_obs,1:3),geom_bndry_panels(ind_obs,4), ...
            geom_bndry_panels(ind_obs,7),geom_bndry_panels(ind_obs,8),...
            geom_bndry_panels(ind_obs,9),smpl_pnts,smpl_wghts,num_diff_pnts);
    end
end

eps0=8.854187817e-12;    
one_over_4pieps0 = (1/(4*pi*eps0));
G_block=G_block * one_over_4pieps0;