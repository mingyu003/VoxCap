function [fN_ch_all]=precond_circulant(dx,L,M,N,Eps_inout)

fl_diel=0;
if (sum(abs(Eps_inout(:,1))) > 1e-13) % check whether we have a dielectric panel
    fl_diel = 1; % yes, we have
end

% constants
eps0=8.854187817e-12;
one_over_4pieps0 = (1/(4*pi*eps0));
% Numerical integration parameters
% Quadrature will be used for the far-interactions (for the interactions
% between panels 100dx away from eachother)in the hybrid calling
num_quad_pnts = 3; % set to 7 if wanted unnecessarily more accuracy
num_diff_pnts = 6; % set to 4 if wanted unnecessarily more accuracy
[smpl_pnts,smpl_wghts] = weights_points(num_quad_pnts,4,[-1 1;-1 1;-1 1;-1 1]);
smpl_wghts = smpl_wghts*(dx^4)*(1/16); smpl_pnts = smpl_pnts*(dx)*(1/2);
len=length(smpl_wghts);

if (fl_diel == 1) % conductor + diel case
    fN_ch_all=cell(3,3,2);
else % only conductor case
    fN_ch_all=cell(3,3,1);
end

dx_vect=[dx dx dx];

% Gxx interactions
tic
G_mn_xx_c = zeros(L+1,M,N);if (fl_diel == 1);G_mn_xx_e = zeros(L+1,M,N);end
n_unit_src=[1 0 0]; n_unit_obs=[1 0 0]; cen_src=[0 0 0];
for mx = 1:L+1
    for my = 1:M
        for mz = 1:N
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_xx_c(mx,my,mz)] = compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len); % for cond
            if (fl_diel == 1)
                if mx==1 && my==1 && mz==1
                    G_mn_xx_e(mx,my,mz)=(dx^2)*2*pi*((Eps_inout(1,1)+Eps_inout(1,2))/(Eps_inout(1,1)-Eps_inout(1,2)));
                else
                    [G_mn_xx_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                        len,num_diff_pnts);
                end
            end
        end
    end
end
G_mn_xx_c = G_mn_xx_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_xx_e = G_mn_xx_e * one_over_4pieps0;end
[fN_ch_all{1,1,1}] = precond_Toeplitz('G_xx_yy_zz_charge',L+1,M,N,G_mn_xx_c); % Gp_mn_xx
clear G_mn_xx_c

if (fl_diel == 1)
    [fN_ch_all{1,1,2}] = precond_Toeplitz('G_xx_Efield',L+1,M,N,G_mn_xx_e); % Gp_mn_xx
    clear G_mn_xx_e
end

disp(['Time for computing circulant tensor of Gxx ::: ',num2str(toc)])
%%%---------------------------------------------------

% Gyy interactions
tic
G_mn_yy_c = zeros(L,M+1,N); if (fl_diel == 1);G_mn_yy_e = zeros(L,M+1,N);end
n_unit_src=[0 1 0]; n_unit_obs=[0 1 0];cen_src=[0 0 0];
for mx = 1:L
    for my = 1:M+1
        for mz = 1:N
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_yy_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
            if (fl_diel == 1)
                if mx==1 && my==1 && mz==1
                    G_mn_yy_e(mx,my,mz)=(dx^2)*2*pi*((Eps_inout(1,1)+Eps_inout(1,2))/(Eps_inout(1,1)-Eps_inout(1,2)));
                else
                    [G_mn_yy_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                        len,num_diff_pnts);
                end
            end
        end
    end
end
G_mn_yy_c = G_mn_yy_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_yy_e = G_mn_yy_e * one_over_4pieps0;end

[fN_ch_all{2,2,1}] = precond_Toeplitz('G_xx_yy_zz_charge',L,M+1,N,G_mn_yy_c); % Gp_mn_yy
clear G_mn_yy_c
if (fl_diel == 1)
    [fN_ch_all{2,2,2}] = precond_Toeplitz('G_yy_Efield',L,M+1,N,G_mn_yy_e); % Gp_mn_yy
    clear G_mn_yy_e
end
disp(['Time for computing circulant tensor of Gyy ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gzz interactions
tic
G_mn_zz_c = zeros(L,M,N+1);if (fl_diel == 1);G_mn_zz_e = zeros(L,M,N+1);end
n_unit_src=[0 0 1]; n_unit_obs=[0 0 1];cen_src=[0 0 0];
for mx = 1:L
    for my = 1:M
        for mz = 1:N+1
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_zz_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
            if (fl_diel == 1)
                if mx==1 && my==1 && mz==1
                    G_mn_zz_e(mx,my,mz)=(dx^2)*2*pi*((Eps_inout(1,1)+Eps_inout(1,2))/(Eps_inout(1,1)-Eps_inout(1,2)));
                else
%                     [G_mn_zz_e(mx,my,mz)] = compute_1over_R_partial_differentiation_integral(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
%                         num_diff_pnts,Eps_inout(1,1),Eps_inout(1,2));
                    [G_mn_zz_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                        len,num_diff_pnts);
                end
            end
        end
    end
end
G_mn_zz_c = G_mn_zz_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_zz_e = G_mn_zz_e * one_over_4pieps0;end
[fN_ch_all{3,3,1}] = precond_Toeplitz('G_xx_yy_zz_charge',L,M,N+1,G_mn_zz_c); % Gp_mn_zz
clear G_mn_zz_c
if (fl_diel == 1)
    [fN_ch_all{3,3,2}] = precond_Toeplitz('G_zz_Efield',L,M,N+1,G_mn_zz_e); % Gp_mn_zz
    clear G_mn_zz_e
end

disp(['Time for computing circulant tensor of Gzz ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gxy interactions
tic
G_mn_xy_c = zeros(L+2,M+2,N); if (fl_diel == 1); G_mn_xy_e = zeros(L+2,M+2,N); end
n_unit_src=[0 1 0]; n_unit_obs=[1 0 0]; cen_src=[dx/2 -dx/2 0];
for mx = 1:L+2
    for my = 1:M+2
        for mz = 1:N
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_xy_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
            if (fl_diel == 1)
                [G_mn_xy_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
end
G_mn_xy_c = G_mn_xy_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_xy_e = G_mn_xy_e * one_over_4pieps0;end

[fN_ch_all{1,2,1}] = precond_Toeplitz('G_xy_charge',L+1,M+1,N,G_mn_xy_c); % Gp_mn_xy
clear G_mn_xy_c
if (fl_diel == 1)   
    [fN_ch_all{1,2,2}] = precond_Toeplitz('G_xy_Efield',L+1,M+1,N,G_mn_xy_e); % Gp_mn_xy
    clear G_mn_xy_e
end

disp(['Time for computing circulant tensor of Gxy ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gxz interactions

tic
G_mn_xz_c = zeros(L+2,M,N+2);
n_unit_src=[0 0 1]; n_unit_obs=[1 0 0]; cen_src=[dx/2 0 -dx/2];
for mx = 1:L+2
    for my = 1:M
        for mz = 1:N+2
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_xz_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
        end
    end
end

G_mn_xz_c = G_mn_xz_c * one_over_4pieps0;

[fN_ch_all{1,3,1}] = precond_Toeplitz('G_xz_charge',L+1,M,N+1,G_mn_xz_c); % Gp_mn_xz
clear G_mn_xz_c
if (fl_diel == 1) %%%%%%here can be optimized, but time display in 1.2.2 needs to be delete.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! remove this part to 1.2.2
    G_mn_xz_e = zeros(L+2,M,N+2);
    n_unit_src=[0 0 1]; n_unit_obs=[1 0 0]; cen_src=[dx/2 0 -dx/2];
    for mx = 1:L+2
        for my = 1:M
            for mz = 1:N+2
                m = [mx my mz];
                cen_obs=((m-1).*dx_vect);
                [G_mn_xz_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
    G_mn_xz_e = G_mn_xz_e * one_over_4pieps0;
    [fN_ch_all{1,3,2}] = precond_Toeplitz('G_xz_Efield',L+1,M,N+1,G_mn_xz_e); % Gp_mn_xz
    clear G_mn_xz_e
end

disp(['Time for computing circulant tensor of Gxz ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gyz interactions
tic
G_mn_yz_c = zeros(L,M+2,N+2); if (fl_diel == 1); G_mn_yz_e = zeros(L,M+2,N+2); end
n_unit_src=[0 0 1]; n_unit_obs=[0 1 0]; cen_src=[0 dx/2 -dx/2];
for mx = 1:L
    for my = 1:M+2
        for mz = 1:N+2
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_yz_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
            if (fl_diel == 1)
                [G_mn_yz_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
end

G_mn_yz_c = G_mn_yz_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_yz_e = G_mn_yz_e * one_over_4pieps0;end
[fN_ch_all{2,3,1}] = precond_Toeplitz('G_yz_charge',L,M+1,N+1,G_mn_yz_c); % Gp_mn_yz
clear G_mn_yz_c
if (fl_diel == 1)  
    [fN_ch_all{2,3,2}] = precond_Toeplitz('G_yz_Efield',L,M+1,N+1,G_mn_yz_e); % Gp_mn_yz
    clear G_mn_yz_e
end

disp(['Time for computing circulant tensor of Gyz ::: ',num2str(toc)])
%%%---------------------------------------------------
% new ones (involving off-diagonal blocks) here just for dielectric !!!
G_mn_yx_c = zeros(L+2,M+2,N);
n_unit_src=[1 0 0]; n_unit_obs=[0 1 0]; cen_src=[-dx/2 dx/2 0];
for mx = 1:L+2
    for my = 1:M+2
        for mz = 1:N
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_yx_c(mx,my,mz)] =  compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
        end
    end
end
G_mn_yx_c = G_mn_yx_c * one_over_4pieps0;
[fN_ch_all{2,1,1}] = precond_Toeplitz('G_yx_charge',L+1,M+1,N,G_mn_yx_c); % Gp_mn_xy
clear G_mn_yx_c

G_mn_zx_c = zeros(L+2,M,N+2);
n_unit_src=[1 0 0]; n_unit_obs=[0 0 1]; cen_src=[-dx/2 0 dx/2];
for mx = 1:L+2
    for my = 1:M
        for mz = 1:N+2
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_zx_c(mx,my,mz)] = compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
        end
    end
end

G_mn_zx_c = G_mn_zx_c * one_over_4pieps0;


[fN_ch_all{3,1,1}] = precond_Toeplitz('G_zx_charge',L+1,M,N+1,G_mn_zx_c); % Gp_mn_xz
clear G_mn_zx_c
G_mn_zy_c = zeros(L,M+2,N+2);
n_unit_src=[0 1 0]; n_unit_obs=[0 0 1]; cen_src=[0 -dx/2 dx/2];
for mx = 1:L
    for my = 1:M+2
        for mz = 1:N+2
            m = [mx my mz];
            cen_obs=((m-1).*dx_vect);
            [G_mn_zy_c(mx,my,mz)] = compute_1overR_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,len);
        end
    end
end

G_mn_zy_c = G_mn_zy_c * one_over_4pieps0;
[fN_ch_all{3,2,1}] = precond_Toeplitz('G_zy_charge',L,M+1,N+1,G_mn_zy_c); % Gp_mn_yz
clear G_mn_zy_c


if (fl_diel == 1)
    % Gyx interactions
    tic
    
    G_mn_yx_e = zeros(L+2,M+2,N);
    n_unit_src=[1 0 0]; n_unit_obs=[0 1 0]; cen_src=[-dx/2 dx/2 0];
    for mx = 1:L+2
        for my = 1:M+2
            for mz = 1:N
                m = [mx my mz];
                cen_obs=((m-1).*dx_vect);
                [G_mn_yx_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
    G_mn_yx_e = G_mn_yx_e * one_over_4pieps0;
    [fN_ch_all{2,1,2}] = precond_Toeplitz('G_yx_Efield',L+1,M+1,N,G_mn_yx_e); % Gp_mn_xy
    clear G_mn_yx_e
    disp(['Time for computing circulant tensor of Gyx ::: ',num2str(toc)])
    %%%---------------------------------------------------
    % Gzx interactions    
    tic
    G_mn_zx_e = zeros(L+2,M,N+2);
    n_unit_src=[1 0 0]; n_unit_obs=[0 0 1]; cen_src=[-dx/2 0 dx/2];
    for mx = 1:L+2
        for my = 1:M
            for mz = 1:N+2
                m = [mx my mz];
                cen_obs=((m-1).*dx_vect);
                [G_mn_zx_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
    
    G_mn_zx_e = G_mn_zx_e * one_over_4pieps0;
    [fN_ch_all{3,1,2}] = precond_Toeplitz('G_zx_Efield',L+1,M,N+1,G_mn_zx_e); % Gp_mn_xz
    clear G_mn_zx_e
    
    disp(['Time for computing circulant tensor of Gzx ::: ',num2str(toc)])
    %%%---------------------------------------------------
    % Gzy interactions
    tic
    
    G_mn_zy_e = zeros(L,M+2,N+2);
    n_unit_src=[0 1 0]; n_unit_obs=[0 0 1]; cen_src=[0 -dx/2 dx/2];
    for mx = 1:L
        for my = 1:M+2
            for mz = 1:N+2
                m = [mx my mz];
                cen_obs=((m-1).*dx_vect);
                [G_mn_zy_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                    len,num_diff_pnts);
            end
        end
    end
    
    G_mn_zy_e = G_mn_zy_e * one_over_4pieps0;
    [fN_ch_all{3,2,2}] = precond_Toeplitz('G_zy_Efield',L,M+1,N+1,G_mn_zy_e); % Gp_mn_yz
    clear G_mn_zy_e

    disp(['Time for computing circulant tensor of Gzy ::: ',num2str(toc)])
end
fN_ch_all{1,1,1}(L+2,:,:)=[];fN_ch_all{1,1,1}(:,M+1,:)=[];fN_ch_all{1,1,1}(:,:,N+1)=[];
fN_ch_all{2,2,1}(L+1,:,:)=[];fN_ch_all{2,2,1}(:,M+2,:)=[];fN_ch_all{2,2,1}(:,:,N+1)=[];
fN_ch_all{3,3,1}(L+1,:,:)=[];fN_ch_all{3,3,1}(:,M+1,:)=[];fN_ch_all{3,3,1}(:,:,N+2)=[];
fN_ch_all{1,2,1}(L+2,:,:)=[];fN_ch_all{1,2,1}(:,M+2,:)=[];fN_ch_all{1,2,1}(:,:,N+1)=[];
fN_ch_all{2,1,1}(L+2,:,:)=[];fN_ch_all{2,1,1}(:,M+2,:)=[];fN_ch_all{2,1,1}(:,:,N+1)=[];
fN_ch_all{1,3,1}(L+2,:,:)=[];fN_ch_all{1,3,1}(:,M+1,:)=[];fN_ch_all{1,3,1}(:,:,N+2)=[];
fN_ch_all{3,1,1}(L+2,:,:)=[];fN_ch_all{3,1,1}(:,M+1,:)=[];fN_ch_all{3,1,1}(:,:,N+2)=[];
fN_ch_all{2,3,1}(L+1,:,:)=[];fN_ch_all{2,3,1}(:,M+2,:)=[];fN_ch_all{2,3,1}(:,:,N+2)=[];
fN_ch_all{3,2,1}(L+1,:,:)=[];fN_ch_all{3,2,1}(:,M+2,:)=[];fN_ch_all{3,2,1}(:,:,N+2)=[];
if (fl_diel == 1)
    fN_ch_all{1,1,2}(L+2,:,:)=[];fN_ch_all{1,1,2}(:,M+1,:)=[];fN_ch_all{1,1,2}(:,:,N+1)=[];
    fN_ch_all{2,2,2}(L+1,:,:)=[];fN_ch_all{2,2,2}(:,M+2,:)=[];fN_ch_all{2,2,2}(:,:,N+1)=[];
    fN_ch_all{3,3,2}(L+1,:,:)=[];fN_ch_all{3,3,2}(:,M+1,:)=[];fN_ch_all{3,3,2}(:,:,N+2)=[];
    fN_ch_all{1,2,2}(L+2,:,:)=[];fN_ch_all{1,2,2}(:,M+2,:)=[];fN_ch_all{1,2,2}(:,:,N+1)=[];
    fN_ch_all{2,1,2}(L+2,:,:)=[];fN_ch_all{2,1,2}(:,M+2,:)=[];fN_ch_all{2,1,2}(:,:,N+1)=[];
    fN_ch_all{1,3,2}(L+2,:,:)=[];fN_ch_all{1,3,2}(:,M+1,:)=[];fN_ch_all{1,3,2}(:,:,N+2)=[];
    fN_ch_all{3,1,2}(L+2,:,:)=[];fN_ch_all{3,1,2}(:,M+1,:)=[];fN_ch_all{3,1,2}(:,:,N+2)=[];
    fN_ch_all{2,3,2}(L+1,:,:)=[];fN_ch_all{2,3,2}(:,M+2,:)=[];fN_ch_all{2,3,2}(:,:,N+2)=[];
    fN_ch_all{3,2,2}(L+1,:,:)=[];fN_ch_all{3,2,2}(:,M+2,:)=[];fN_ch_all{3,2,2}(:,:,N+2)=[];
end
% %----------------------------------------------------------------------------------

disp('-----------------------------------------------------')