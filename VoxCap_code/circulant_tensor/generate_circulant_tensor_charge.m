function [fN_charge]=generate_circulant_tensor_charge(dx,L,M,N,Eps_inout,tolerance,fl_Tucker_decomp)
% This routine computes the circulant tensors for fast computation of
% potentials due to charges on square panels. Two different computation
% schemes are implemented: (i) circulant from Toeplitz and (ii) circulant
% via direct computation
%%%INPUTS
%dx: resolution of voxels
%L,M,N: number of voxels in x,y and z direction, respectively, in
%       computational domain
%tolerance: tolerancefor Turker compression
%Eps_inout: [eps_in, eps_out]
fl_circ_from_Toep = 1; % 1 for computing circulant from Toeplitz embedding
fl_no_fft = 0; % 1 if no FFT is applied to circulant tensors

disp('-----------------------------------------------------')
if (fl_circ_from_Toep == 1)
    disp('Computing circulant tensors via embedding their Toeplitz')
else
    disp('Computing circulant tensors (all entries) ')
end

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

if (fl_diel == 1) % conductor + diel case
    
    num_elem_est=(2*(L+1)*(M+1)*(N+1)); memestimated=num_elem_est*8/(1024*1024);
    disp(['Memory for temporarily storing a Toeplitz tensor (MB) ::: ' , num2str(memestimated)]);
    
    num_elem_est=15*(2*(L+1)*(2*(M+1))*(2*(N+1))); memestimated=num_elem_est*16/(1024*1024);
    disp(['Memory for storing circulant tensors (MB) ::: ' , num2str(memestimated)]);
    
else % only conductor case
    
    num_elem_est=((L+1)*(M+1)*(N+1)); memestimated=num_elem_est*8/(1024*1024);
    disp(['Memory for temporarily storing a Toeplitz tensor (MB) ::: ' , num2str(memestimated)]);
    
    num_elem_est=6*(2*(L+1)*(2*(M+1))*(2*(N+1))); memestimated=num_elem_est*16/(1024*1024);
    disp(['Memory for storing circulant tensors (MB) ::: ' , num2str(memestimated)]);    
end

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
                    G_mn_xx_e(mx,my,mz)=0;
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
[fN_ch_all{1,1,1}] = Toeplitz_JVIE('G_xx_charge',L+1,M+1,N+1,G_mn_xx_c); % Gp_mn_xx
clear G_mn_xx_c

if (fl_diel == 1)   
    [fN_ch_all{1,1,2}] = Toeplitz_JVIE('G_xx_Efield',L+1,M+1,N+1,G_mn_xx_e); % Gp_mn_xx
    clear G_mn_xx_e
end

disp(['Time for computing circulant tensor of Gxx ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if fl_Tucker_decomp==1
    if (fl_diel == 1)
        fN_charge=cell(73,1);
    else
        fN_charge=cell(37,1);
    end
end
if (fl_no_fft == 0)
    tic
    fN_ch_all{1,1,1}(:,:,:) = fftn(fN_ch_all{1,1,1}(:,:,:));
    if (fl_diel == 1)
        fN_ch_all{1,1,2}(:,:,:) = fftn(fN_ch_all{1,1,2}(:,:,:));
    end
end
disp(['Time for FFT of Gxx ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
    [c_xx_factor_matrix1,c_xx_factor_matrix2,c_xx_factor_matrix3,c_xx_core_tensor,c_xx_mem] = Tucker(fN_ch_all{1,1,1},tolerance);
    fN_ch_all{1,1,1}=[];
    fN_charge{2}=c_xx_core_tensor;
    fN_charge{3}=c_xx_factor_matrix1;
    fN_charge{4}=c_xx_factor_matrix2;
    fN_charge{5}=c_xx_factor_matrix3;% condtctor_xx
    clear c_xx_factor_matrix1 c_xx_factor_matrix2 c_xx_factor_matrix3 c_xx_core_tensor   ;
    if (fl_diel == 1)
        [d_xx_factor_matrix1,d_xx_factor_matrix2,d_xx_factor_matrix3,d_xx_core_tensor,d_xx_mem] = Tucker(fN_ch_all{1,1,2},tolerance);
        fN_ch_all{1,1,2}=[];
        fN_charge{38}=d_xx_core_tensor;% dielectric_xx
        fN_charge{39}=d_xx_factor_matrix1;
        fN_charge{40}=d_xx_factor_matrix2;
        fN_charge{41}=d_xx_factor_matrix3;
        clear d_xx_factor_matrix1 d_xx_factor_matrix2 d_xx_factor_matrix3 d_xx_core_tensor  ;
    end
    disp(['Time for compression of Gxx ::: ',num2str(toc)])
end

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
                    G_mn_yy_e(mx,my,mz)=0;
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
[fN_ch_all{2,2,1}] = Toeplitz_JVIE('G_yy_charge',L+1,M+1,N+1,G_mn_yy_c); % Gp_mn_yy
clear G_mn_yy_c
if (fl_diel == 1)
    [fN_ch_all{2,2,2}] = Toeplitz_JVIE('G_yy_Efield',L+1,M+1,N+1,G_mn_yy_e); % Gp_mn_yy
    clear G_mn_yy_e
end
disp(['Time for computing circulant tensor of Gyy ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if (fl_no_fft == 0)
    tic
    fN_ch_all{2,2,1}(:,:,:) = fftn(fN_ch_all{2,2,1}(:,:,:));
    if (fl_diel == 1)
        fN_ch_all{2,2,2}(:,:,:) = fftn(fN_ch_all{2,2,2}(:,:,:));
    end
end
disp(['Time for FFT of Gyy ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_y, MfN_y, NfN_y, ~] = size(fN_ch_all{2,2,1}); %yy
    [c_yy_factor_matrix1,c_yy_factor_matrix2,c_yy_factor_matrix3,c_yy_core_tensor,c_yy_mem] = Tucker(fN_ch_all{2,2,1},tolerance);
    fN_ch_all{2,2,1}=[];
    fN_charge{6}=c_yy_core_tensor;
    fN_charge{7}=c_yy_factor_matrix1;
    fN_charge{8}=c_yy_factor_matrix2;
    fN_charge{9}=c_yy_factor_matrix3;% condtctor_yy
    clear c_yy_factor_matrix1 c_yy_factor_matrix2 c_yy_factor_matrix3 c_yy_core_tensor   ;
    if (fl_diel == 1)
        [d_yy_factor_matrix1,d_yy_factor_matrix2,d_yy_factor_matrix3,d_yy_core_tensor,d_yy_mem] = Tucker(fN_ch_all{2,2,2},tolerance);
        fN_ch_all{2,2,2}=[];
        fN_charge{42}=d_yy_core_tensor;% dielectric_yy
        fN_charge{43}=d_yy_factor_matrix1;
        fN_charge{44}=d_yy_factor_matrix2;
        fN_charge{45}=d_yy_factor_matrix3;
        clear d_yy_factor_matrix1 d_yy_factor_matrix2 d_yy_factor_matrix3 d_yy_core_tensor  ;
    end
    disp(['Time for compression of Gyy ::: ',num2str(toc)])
end

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
                    G_mn_zz_e(mx,my,mz)=0;
                else
                    [G_mn_zz_e(mx,my,mz)] = compute_1overR_partial_diff_integral_mexn(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,...
                        len,num_diff_pnts);
%                     [G_mn_zz_e(mx,my,mz)] =compute_1over_R_partial_differentiation_integral(dx,cen_src,cen_obs,n_unit_src,n_unit_obs,smpl_pnts,smpl_wghts,num_diff_pnts,2,1);

                end
            end
        end
    end
end
G_mn_zz_c = G_mn_zz_c * one_over_4pieps0;
if (fl_diel == 1);G_mn_zz_e = G_mn_zz_e * one_over_4pieps0;end
[fN_ch_all{3,3,1}] = Toeplitz_JVIE('G_zz_charge',L+1,M+1,N+1,G_mn_zz_c); % Gp_mn_zz
clear G_mn_zz_c
if (fl_diel == 1)      
    [fN_ch_all{3,3,2}] = Toeplitz_JVIE('G_zz_Efield',L+1,M+1,N+1,G_mn_zz_e); % Gp_mn_zz
    clear G_mn_zz_e
end

disp(['Time for computing circulant tensor of Gzz ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if (fl_no_fft == 0)
    tic
    fN_ch_all{3,3,1}(:,:,:) = fftn(fN_ch_all{3,3,1}(:,:,:));
    if (fl_diel == 1)
        fN_ch_all{3,3,2}(:,:,:) = fftn(fN_ch_all{3,3,2}(:,:,:));
    end
end
disp(['Time for FFT of Gzz ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_z, MfN_z, NfN_z, ~] = size(fN_ch_all{3,3,1}); %zz
    [c_zz_factor_matrix1,c_zz_factor_matrix2,c_zz_factor_matrix3,c_zz_core_tensor,c_zz_mem] = Tucker(fN_ch_all{3,3,1},tolerance);
    fN_ch_all{3,3,1}=[];
    fN_charge{10}=c_zz_core_tensor;
    fN_charge{11}=c_zz_factor_matrix1;
    fN_charge{12}=c_zz_factor_matrix2;
    fN_charge{13}=c_zz_factor_matrix3;% condtctor_zz
    clear c_zz_factor_matrix1 c_zz_factor_matrix2 c_zz_factor_matrix3 c_zz_core_tensor   ;
    if (fl_diel == 1)
        [d_zz_factor_matrix1,d_zz_factor_matrix2,d_zz_factor_matrix3,d_zz_core_tensor,d_zz_mem] = Tucker(fN_ch_all{3,3,2},tolerance);
        fN_ch_all{3,3,2}=[];
        fN_charge{46}=d_zz_core_tensor;% dielectric_zz
        fN_charge{47}=d_zz_factor_matrix1;
        fN_charge{48}=d_zz_factor_matrix2;
        fN_charge{49}=d_zz_factor_matrix3;
        clear d_zz_factor_matrix1 d_zz_factor_matrix2 d_zz_factor_matrix3 d_zz_core_tensor  ;
    end
    disp(['Time for compression of Gzz ::: ',num2str(toc)])
end

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
[fN_ch_all{1,2,1}] = Toeplitz_JVIE('G_xy_charge',L+1,M+1,N+1,G_mn_xy_c); % Gp_mn_xy
clear G_mn_xy_c
if (fl_diel == 1)   
    [fN_ch_all{1,2,2}] = Toeplitz_JVIE('G_xy_Efield',L+1,M+1,N+1,G_mn_xy_e); % Gp_mn_xy
    clear G_mn_xy_e
end

disp(['Time for computing circulant tensor of Gxy ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if (fl_no_fft == 0)
    tic
    fN_ch_all{1,2,1}(:,:,:) = fftn(fN_ch_all{1,2,1}(:,:,:));
    if (fl_diel == 1)
        fN_ch_all{1,2,2}(:,:,:) = fftn(fN_ch_all{1,2,2}(:,:,:));
    end
end
disp(['Time for FFT of Gxy ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_xy, MfN_xy, NfN_xy, ~] = size(fN_ch_all{1,2,1}); %xy
    [c_xy_factor_matrix1,c_xy_factor_matrix2,c_xy_factor_matrix3,c_xy_core_tensor,c_xy_mem] = Tucker(fN_ch_all{1,2,1},tolerance);
    [c_yx_factor_matrix1,c_yx_factor_matrix2,c_yx_factor_matrix3,c_yx_core_tensor,c_yx_mem] = Tucker(conj(fN_ch_all{1,2,1}),tolerance);
    fN_ch_all{1,2,1}=[];
    fN_charge{14}=c_xy_core_tensor;
    fN_charge{15}=c_xy_factor_matrix1;
    fN_charge{16}=c_xy_factor_matrix2;
    fN_charge{17}=c_xy_factor_matrix3;% condtctor_xy
    fN_charge{26}=c_yx_core_tensor;
    fN_charge{27}=c_yx_factor_matrix1;
    fN_charge{28}=c_yx_factor_matrix2;
    fN_charge{29}=c_yx_factor_matrix3;% condtctor_yx
    clear c_xy_factor_matrix1 c_xy_factor_matrix2 c_xy_factor_matrix3 c_xy_core_tensor c_yx_factor_matrix1 c_yx_factor_matrix2 c_yx_factor_matrix3 c_yx_core_tensor;
    
    
    if (fl_diel == 1)
        [d_xy_factor_matrix1,d_xy_factor_matrix2,d_xy_factor_matrix3,d_xy_core_tensor,d_xy_mem] = Tucker(fN_ch_all{1,2,2},tolerance);
%         fN_ch_all{1,2,2}=[];
        fN_charge{50}=d_xy_core_tensor;% dielectric_xy
        fN_charge{51}=d_xy_factor_matrix1;
        fN_charge{52}=d_xy_factor_matrix2;
        fN_charge{53}=d_xy_factor_matrix3;
        clear d_xy_factor_matrix1 d_xy_factor_matrix2 d_xy_factor_matrix3 d_xy_core_tensor  ;
    end
    disp(['Time for compression of Gxy ::: ',num2str(toc)])
end

%%%---------------------------------------------------
% Gxz interactions

tic
G_mn_xz_c = zeros(L+2,M,N+2);% if (fl_diel == 1); G_mn_xz_e = zeros(L+2,M,N+2); end
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
[fN_ch_all{1,3,1}] = Toeplitz_JVIE('G_xz_charge',L+1,M+1,N+1,G_mn_xz_c); % Gp_mn_xz
clear G_mn_xz_c
if (fl_diel == 1) %%%%%%here can be optimized, but time display in 1.2.2 needs to be delete.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! remove this part to 1.2.2
    if N==M
        for ii=1:2*(M+1)
            for jj=1:2*(N+1)
                fN_ch_all{1,3,2}(:,ii,jj)=fN_ch_all{1,2,2}(:,jj,ii);
            end
        end
    else
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
        [fN_ch_all{1,3,2}] = Toeplitz_JVIE('G_xz_Efield',L+1,M+1,N+1,G_mn_xz_e); % Gp_mn_xz
        clear G_mn_xz_e
    end
end

disp(['Time for computing circulant tensor of Gxz ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if (fl_no_fft == 0)
    tic
    fN_ch_all{1,3,1}(:,:,:) = fftn(fN_ch_all{1,3,1}(:,:,:));
    if (fl_diel == 1) && abs(M-N)>1e-12
        fN_ch_all{1,3,2}(:,:,:) = fftn(fN_ch_all{1,3,2}(:,:,:));
    end
end
disp(['Time for FFT of Gxz ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
    [c_xz_factor_matrix1,c_xz_factor_matrix2,c_xz_factor_matrix3,c_xz_core_tensor,c_xz_mem] = Tucker(fN_ch_all{1,3,1},tolerance);
    [c_zx_factor_matrix1,c_zx_factor_matrix2,c_zx_factor_matrix3,c_zx_core_tensor,c_zx_mem] = Tucker(conj(fN_ch_all{1,3,1}),tolerance);
    fN_ch_all{1,3,1}=[];
    fN_charge{18}=c_xz_core_tensor;
    fN_charge{19}=c_xz_factor_matrix1;
    fN_charge{20}=c_xz_factor_matrix2;
    fN_charge{21}=c_xz_factor_matrix3;% condtctor_xz
    fN_charge{30}=c_zx_core_tensor;
    fN_charge{31}=c_zx_factor_matrix1;
    fN_charge{32}=c_zx_factor_matrix2;
    fN_charge{33}=c_zx_factor_matrix3;% condtctor_zx
    clear c_xz_factor_matrix1 c_xz_factor_matrix2 c_xz_factor_matrix3 c_xz_core_tensor c_zx_factor_matrix1 c_zx_factor_matrix2 c_zx_factor_matrix3 c_zx_core_tensor ;
    clear c_xz_factor_matrix1 c_xz_factor_matrix2 c_xz_factor_matrix3 c_xz_core_tensor   ;
    if (fl_diel == 1)
        [d_xz_factor_matrix1,d_xz_factor_matrix2,d_xz_factor_matrix3,d_xz_core_tensor,d_xz_mem] = Tucker(fN_ch_all{1,3,2},tolerance);
%         fN_ch_all{1,3,2}=[];
        fN_charge{54}=d_xz_core_tensor;% dielectric_xz
        fN_charge{55}=d_xz_factor_matrix1;
        fN_charge{56}=d_xz_factor_matrix2;
        fN_charge{57}=d_xz_factor_matrix3;
        clear d_xz_factor_matrix1 d_xz_factor_matrix2 d_xz_factor_matrix3 d_xz_core_tensor  ;
    end
    disp(['Time for compression of Gxz ::: ',num2str(toc)])
end

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
[fN_ch_all{2,3,1}] = Toeplitz_JVIE('G_yz_charge',L+1,M+1,N+1,G_mn_yz_c); % Gp_mn_yz
clear G_mn_yz_c
if (fl_diel == 1)
    [fN_ch_all{2,3,2}] = Toeplitz_JVIE('G_yz_Efield',L+1,M+1,N+1,G_mn_yz_e); % Gp_mn_yz
    clear G_mn_yz_e
end

disp(['Time for computing circulant tensor of Gyz ::: ',num2str(toc)])
%%%----------------------------------------------
%               FFT & Tucker
%%%----------------------------------------------
if (fl_no_fft == 0)
    tic
    fN_ch_all{2,3,1}(:,:,:) = fftn(fN_ch_all{2,3,1}(:,:,:));
    if (fl_diel == 1)
        fN_ch_all{2,3,2}(:,:,:) = fftn(fN_ch_all{2,3,2}(:,:,:));
    end
end
disp(['Time for FFT of Gyz ::: ',num2str(toc)])
if fl_Tucker_decomp==1
    tic
    [LfN_yz, MfN_yz, NfN_yz, ~] = size(fN_ch_all{2,3,1}); %yz
    [c_yz_factor_matrix1,c_yz_factor_matrix2,c_yz_factor_matrix3,c_yz_core_tensor,c_yz_mem] = Tucker(fN_ch_all{2,3,1},tolerance);
    [c_zy_factor_matrix1,c_zy_factor_matrix2,c_zy_factor_matrix3,c_zy_core_tensor,c_zy_mem] = Tucker(conj(fN_ch_all{2,3,1}),tolerance);
    fN_ch_all{2,3,1}=[];
    fN_charge{22}=c_yz_core_tensor;
    fN_charge{23}=c_yz_factor_matrix1;
    fN_charge{24}=c_yz_factor_matrix2;
    fN_charge{25}=c_yz_factor_matrix3;% condtctor_yz
    fN_charge{34}=c_zy_core_tensor;
    fN_charge{35}=c_zy_factor_matrix1;
    fN_charge{36}=c_zy_factor_matrix2;
    fN_charge{37}=c_zy_factor_matrix3;% condtctor_zy
    clear c_yz_factor_matrix1 c_yz_factor_matrix2 c_yz_factor_matrix3 c_yz_core_tensor c_zy_factor_matrix1 c_zy_factor_matrix2 c_zy_factor_matrix3 c_zy_core_tensor;
    if (fl_diel == 1)
        [d_yz_factor_matrix1,d_yz_factor_matrix2,d_yz_factor_matrix3,d_yz_core_tensor,d_yz_mem] = Tucker(fN_ch_all{2,3,2},tolerance);
%         fN_ch_all{2,3,2}=[];
        fN_charge{58}=d_yz_core_tensor;% dielectric_yz
        fN_charge{59}=d_yz_factor_matrix1;
        fN_charge{60}=d_yz_factor_matrix2;
        fN_charge{61}=d_yz_factor_matrix3;
        clear d_yz_factor_matrix1 d_yz_factor_matrix2 d_yz_factor_matrix3 d_yz_core_tensor  ;
    end
    disp(['Time for compression of Gyz ::: ',num2str(toc)])
end

%%%---------------------------------------------------
% new ones (involving off-diagonal blocks) here just for dielectric !!!
if (fl_diel == 1)
    % Gyx interactions
    tic
    if L==M
        for ii=1:2*(L+1)
            for jj=1:2*(M+1)
                fN_ch_all{2,1,2}(ii,jj,:)=fN_ch_all{1,2,2}(jj,ii,:);% maybe there are bugs!!!!!!!
            end
        end
    else
        
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
        [fN_ch_all{2,1,2}] = Toeplitz_JVIE('G_yx_Efield',L+1,M+1,N+1,G_mn_yx_e); % Gp_mn_xy
        clear G_mn_yx_e
    end
    disp(['Time for computing circulant tensor of Gyx ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        if abs(L-M)>1e-12
            fN_ch_all{2,1,2}(:,:,:) = fftn(fN_ch_all{2,1,2}(:,:,:));
        end
    end
    disp(['Time for FFT of Gyx ::: ',num2str(toc)])
    if fl_Tucker_decomp==1
        tic
        [d_yx_factor_matrix1,d_yx_factor_matrix2,d_yx_factor_matrix3,d_yx_core_tensor,d_yx_mem] = Tucker(fN_ch_all{2,1,2},tolerance);
        fN_ch_all{2,1,2}=[];
        fN_charge{62}=d_yx_core_tensor;% dielectric_yx
        fN_charge{63}=d_yx_factor_matrix1;
        fN_charge{64}=d_yx_factor_matrix2;
        fN_charge{65}=d_yx_factor_matrix3;
        clear d_yx_factor_matrix1 d_yx_factor_matrix2 d_yx_factor_matrix3 d_yx_core_tensor  ;
        disp(['Time for compression of Gyx ::: ',num2str(toc)])
    end
    
    %%%---------------------------------------------------
    % Gzx interactions
    tic
    if L==N
        for ii=1:2*(L+1)
            for jj=1:2*(N+1)
                fN_ch_all{3,1,2}(ii,:,jj)=fN_ch_all{1,3,2}(jj,:,ii);% maybe there are bugs!!!!!!!
            end
        end
    else
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
        [fN_ch_all{3,1,2}] = Toeplitz_JVIE('G_zx_Efield',L+1,M+1,N+1,G_mn_zx_e); % Gp_mn_xz
        clear G_mn_zx_e
    end
    disp(['Time for computing circulant tensor of Gzx ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        if abs(L-N)>1e-12
            fN_ch_all{3,1,2}(:,:,:) = fftn(fN_ch_all{3,1,2}(:,:,:));
        end
    end
    disp(['Time for FFT of Gzx ::: ',num2str(toc)])
    if fl_Tucker_decomp==1
        tic
        [d_zx_factor_matrix1,d_zx_factor_matrix2,d_zx_factor_matrix3,d_zx_core_tensor,d_zx_mem] = Tucker(fN_ch_all{3,1,2},tolerance);
        fN_ch_all{3,1,2}=[];
        fN_charge{66}=d_zx_core_tensor;% dielectric_zx
        fN_charge{67}=d_zx_factor_matrix1;
        fN_charge{68}=d_zx_factor_matrix2;
        fN_charge{69}=d_zx_factor_matrix3;
        clear d_zx_factor_matrix1 d_zx_factor_matrix2 d_zx_factor_matrix3 d_zx_core_tensor  ;
        disp(['Time for compression of Gzx ::: ',num2str(toc)])
    end
    
    %%%---------------------------------------------------
    % Gzy interactions
    tic
    if M==N
        for ii=1:2*(M+1)
            for jj=1:2*(N+1)
                fN_ch_all{3,2,2}(:,ii,jj)=fN_ch_all{2,3,2}(:,jj,ii);% maybe there are bugs!!!!!!!
            end
        end
    else
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
        [fN_ch_all{3,2,2}] = Toeplitz_JVIE('G_zy_Efield',L+1,M+1,N+1,G_mn_zy_e); % Gp_mn_yz
        clear G_mn_zy_e
    end
    
    disp(['Time for computing circulant tensor of Gzy ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        if abs(M-N)>1e-12
            fN_ch_all{3,2,2}(:,:,:) = fftn(fN_ch_all{3,2,2}(:,:,:));
        end
    end
    disp(['Time for FFT of Gzy ::: ',num2str(toc)])
    if fl_Tucker_decomp==1
        tic
        [d_zy_factor_matrix1,d_zy_factor_matrix2,d_zy_factor_matrix3,d_zy_core_tensor,d_zy_mem] = Tucker(fN_ch_all{3,2,2},tolerance);
        fN_ch_all{3,2,2}=[];
        fN_charge{70}=d_zy_core_tensor;% dielectric_zy
        fN_charge{71}=d_zy_factor_matrix1;
        fN_charge{72}=d_zy_factor_matrix2;
        fN_charge{73}=d_zy_factor_matrix3;
        clear d_zy_factor_matrix1 d_zy_factor_matrix2 d_zy_factor_matrix3 d_zy_core_tensor  ;
        disp(['Time for compression of Gzy ::: ',num2str(toc)])
    end
    
    %%%---------------------------------------------------
end
if fl_Tucker_decomp==1
    fN_size=[LfN_x, MfN_x, NfN_x;LfN_y, MfN_y, NfN_y;LfN_z, MfN_z, NfN_z;LfN_xy, MfN_xy, NfN_xy;LfN_xz, MfN_xz, NfN_xz;LfN_yz, MfN_yz, NfN_yz];
    fN_charge{1}=fN_size;%fN_ch_all size
end

if fl_Tucker_decomp==1
    
   if (fl_diel == 1)
        memory_compress_diel=(c_xx_mem+c_yy_mem+c_zz_mem+c_xy_mem+c_yx_mem+c_xz_mem+c_zx_mem+c_yz_mem+c_zy_mem+d_xx_mem+d_yy_mem+d_zz_mem+d_xy_mem+d_xz_mem+d_yz_mem+d_yx_mem+d_zx_mem+d_zy_mem)*16/1024^2;
        disp(['Memory for storing cross_Tucker-compressed circulant tensors (MB) ::: ' , num2str(memory_compress_diel)]);
    else
        memory_compress=(c_xx_mem+c_yy_mem+c_zz_mem+c_xy_mem+c_yx_mem+c_xz_mem+c_zx_mem+c_yz_mem+c_zy_mem)*16/1024^2;
        disp(['Memory for storing cross_Tucker-compressed circulant tensors (MB) ::: ' , num2str(memory_compress)]);
    end
else
    fN_charge=fN_ch_all;
end
%-------------------------------------------------------------------------
    
% %----------------------------------------------------------------------------------
if (fl_circ_from_Toep == 1)
    disp('Done... Computing circulant tensors via embedding their Toeplitz')
else
    disp('Done... Computing circulant tensors (all entries) ')
end
disp('-----------------------------------------------------')