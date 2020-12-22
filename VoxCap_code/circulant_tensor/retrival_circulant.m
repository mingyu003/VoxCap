function [fN_charge]=retrival_circulant(dx,L,M,N,Eps_inout,tolerance,fl_Tucker_decomp,fN_charge_Toeplitz)
% This routine computes the circulant tensors by retrieval prestored
% non-FFT'ed Toeplitz tensor.

%%%INPUTS
%dx: resolution of voxels
%L,M,N: number of voxels in x,y and z direction, respectively, in
%       computational domain
%tolerance: tolerancefor Turker compression
%Eps_inout: [eps_in, eps_out]
%fN_charge_Toeplitz: prestored non-FFT'ed Toeplitz tensor

size_stored_tensor=size(fN_charge_Toeplitz);
fl_no_fft = 0; % 1 if no FFT is applied to circulant tensors
disp('-----------------------------------------------------')
disp('Computing circulant tensors via embedding their Toeplitz')
fl_diel=0;
if (sum(abs(Eps_inout(:,1))) > 1e-13) % check whether we have a dielectric panel
    fl_diel = 1; % yes, we have
end

if (fl_diel == 1) % conductor + diel case
    fN_ch_all=cell(3,3,2);
else % only conductor case
    fN_ch_all=cell(3,3,1);
end

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
if fl_Tucker_decomp==1
    if (fl_diel == 1)
        fN_charge=cell(73,1);
    else
        fN_charge=cell(37,1);
    end
end
% Gxx interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_xx_c = ten_mat_prod(fN_charge_Toeplitz{2},{fN_charge_Toeplitz{3},fN_charge_Toeplitz{4},fN_charge_Toeplitz{5}});
    [fN_ch_all{1,1,1}] = Toeplitz_JVIE('G_xx_charge',L+1,M+1,N+1,dx^3*G_mn_xx_c(1:L+1,1:M,1:N)); % Gp_mn_xx
    clear G_mn_xx_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,1,1}(:,:,:) = fftn(fN_ch_all{1,1,1}(:,:,:));
        disp(['Time for FFT of Gxx_c ::: ',num2str(toc)])
        if fl_Tucker_decomp==1
            tic
            [LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
            [c_xx_factor_matrix1,c_xx_factor_matrix2,c_xx_factor_matrix3,c_xx_core_tensor,c_xx_mem] = Tucker(fN_ch_all{1,1,1},tolerance);
            fN_ch_all{1,1,1}=[];
            fN_charge{2}=c_xx_core_tensor;
            fN_charge{3}=c_xx_factor_matrix1;
            fN_charge{4}=c_xx_factor_matrix2;
            fN_charge{5}=c_xx_factor_matrix3;% condtctor_xx
            clear c_xx_factor_matrix1 c_xx_factor_matrix2 c_xx_factor_matrix3 c_xx_core_tensor;
            disp(['Time for compression of Gxx_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_xx_e = ten_mat_prod(fN_charge_Toeplitz{26},{fN_charge_Toeplitz{27},fN_charge_Toeplitz{28},fN_charge_Toeplitz{29}});
        [fN_ch_all{1,1,2}] = Toeplitz_JVIE('G_xx_Efield',L+1,M+1,N+1,dx^2*G_mn_xx_e(1:L+1,1:M,1:N)); % Gp_mn_xx
        clear G_mn_xx_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{1,1,2}(:,:,:) = fftn(fN_ch_all{1,1,2}(:,:,:));
            disp(['Time for FFT of Gxx_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xx_factor_matrix1,d_xx_factor_matrix2,d_xx_factor_matrix3,d_xx_core_tensor,d_xx_mem] = Tucker(fN_ch_all{1,1,2},tolerance);
                fN_ch_all{1,1,2}=[];
                fN_charge{38}=d_xx_core_tensor;% dielectric_xx
                fN_charge{39}=d_xx_factor_matrix1;
                fN_charge{40}=d_xx_factor_matrix2;
                fN_charge{41}=d_xx_factor_matrix3;
                clear d_xx_factor_matrix1 d_xx_factor_matrix2 d_xx_factor_matrix3 d_xx_core_tensor  ;
                disp(['Time for compression of Gxx_d ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{1,1,1}] = Toeplitz_JVIE('G_xx_charge',L+1,M+1,N+1,dx^3*fN_charge_Toeplitz{1,1,1}(1:L+1,1:M,1:N)); % Gp_mn_xx
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,1,1}(:,:,:) = fftn(fN_ch_all{1,1,1}(:,:,:));
        disp(['Time for FFT of Gxx_c ::: ',num2str(toc)])
        if fl_Tucker_decomp==1
            tic
            [LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
            [c_xx_factor_matrix1,c_xx_factor_matrix2,c_xx_factor_matrix3,c_xx_core_tensor,c_xx_mem] = Tucker(fN_ch_all{1,1,1},tolerance);
            fN_ch_all{1,1,1}=[];
            fN_charge{2}=c_xx_core_tensor;
            fN_charge{3}=c_xx_factor_matrix1;
            fN_charge{4}=c_xx_factor_matrix2;
            fN_charge{5}=c_xx_factor_matrix3;% condtctor_xx
            clear c_xx_factor_matrix1 c_xx_factor_matrix2 c_xx_factor_matrix3 c_xx_core_tensor;
            disp(['Time for compression of Gxx_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        [fN_ch_all{1,1,2}] = Toeplitz_JVIE('G_xx_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{1,1,2}(1:L+1,1:M,1:N)); % Gp_mn_xx
        if (fl_no_fft == 0)
            tic
            fN_ch_all{1,1,2}(:,:,:) = fftn(fN_ch_all{1,1,2}(:,:,:));
            disp(['Time for FFT of Gxx_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xx_factor_matrix1,d_xx_factor_matrix2,d_xx_factor_matrix3,d_xx_core_tensor,d_xx_mem] = Tucker(fN_ch_all{1,1,2},tolerance);
                fN_ch_all{1,1,2}=[];
                fN_charge{38}=d_xx_core_tensor;% dielectric_xx
                fN_charge{39}=d_xx_factor_matrix1;
                fN_charge{40}=d_xx_factor_matrix2;
                fN_charge{41}=d_xx_factor_matrix3;
                clear d_xx_factor_matrix1 d_xx_factor_matrix2 d_xx_factor_matrix3 d_xx_core_tensor  ;
                disp(['Time for compression of Gxx_d ::: ',num2str(toc)])
            end
        end
    end
end
disp(['Total time for computing circulant tensor of Gxx ::: ',num2str(toc)])

%%%---------------------------------------------------

% Gyy interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_yy_c = ten_mat_prod(fN_charge_Toeplitz{6},{fN_charge_Toeplitz{7},fN_charge_Toeplitz{8},fN_charge_Toeplitz{9}});
    [fN_ch_all{2,2,1}] = Toeplitz_JVIE('G_yy_charge',L+1,M+1,N+1,dx^3*G_mn_yy_c(1:L,1:M+1,1:N)); % Gp_mn_yy
    clear G_mn_yy_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{2,2,1}(:,:,:) = fftn(fN_ch_all{2,2,1}(:,:,:));
        disp(['Time for FFT of Gyy_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gyy_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_yy_e = ten_mat_prod(fN_charge_Toeplitz{30},{fN_charge_Toeplitz{31},fN_charge_Toeplitz{32},fN_charge_Toeplitz{33}});
        [fN_ch_all{2,2,2}] = Toeplitz_JVIE('G_yy_Efield',L+1,M+1,N+1,dx^2*G_mn_yy_e(1:L,1:M+1,1:N)); % Gp_mn_yy
        clear G_mn_yy_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{2,2,2}(:,:,:) = fftn(fN_ch_all{2,2,2}(:,:,:));
            disp(['Time for FFT of Gyy_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_yy_factor_matrix1,d_yy_factor_matrix2,d_yy_factor_matrix3,d_yy_core_tensor,d_yy_mem] = Tucker(fN_ch_all{2,2,2},tolerance);
                fN_ch_all{2,2,2}=[];
                fN_charge{42}=d_yy_core_tensor;% dielectric_yy
                fN_charge{43}=d_yy_factor_matrix1;
                fN_charge{44}=d_yy_factor_matrix2;
                fN_charge{45}=d_yy_factor_matrix3;
                clear d_yy_factor_matrix1 d_yy_factor_matrix2 d_yy_factor_matrix3 d_yy_core_tensor  ;
                disp(['Time for compression of Gyy_d ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{2,2,1}] = Toeplitz_JVIE('G_yy_charge',L+1,M+1,N+1,dx^3.*fN_charge_Toeplitz{2,2,1}(1:L,1:M+1,1:N)); % Gp_mn_yy
    if (fl_no_fft == 0)
        tic
        fN_ch_all{2,2,1}(:,:,:) = fftn(fN_ch_all{2,2,1}(:,:,:));
        disp(['Time for FFT of Gyy_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gyy_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        [fN_ch_all{2,2,2}] = Toeplitz_JVIE('G_yy_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{2,2,2}(1:L,1:M+1,1:N)); % Gp_mn_yy
        if (fl_no_fft == 0)
            tic
            fN_ch_all{2,2,2}(:,:,:) = fftn(fN_ch_all{2,2,2}(:,:,:));
            disp(['Time for FFT of Gyy_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_yy_factor_matrix1,d_yy_factor_matrix2,d_yy_factor_matrix3,d_yy_core_tensor,d_yy_mem] = Tucker(fN_ch_all{2,2,2},tolerance);
                fN_ch_all{2,2,2}=[];
                fN_charge{42}=d_yy_core_tensor;% dielectric_yy
                fN_charge{43}=d_yy_factor_matrix1;
                fN_charge{44}=d_yy_factor_matrix2;
                fN_charge{45}=d_yy_factor_matrix3;
                clear d_yy_factor_matrix1 d_yy_factor_matrix2 d_yy_factor_matrix3 d_yy_core_tensor  ;
                disp(['Time for compression of Gyy_d ::: ',num2str(toc)])
            end
        end
    end
end

disp(['Total time for computing circulant tensor of Gyy ::: ',num2str(toc)])

%%%---------------------------------------------------
% Gzz interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_zz_c = ten_mat_prod(fN_charge_Toeplitz{10},{fN_charge_Toeplitz{11},fN_charge_Toeplitz{12},fN_charge_Toeplitz{13}});
    [fN_ch_all{3,3,1}] = Toeplitz_JVIE('G_zz_charge',L+1,M+1,N+1,dx^3*G_mn_zz_c(1:L,1:M,1:N+1)); % Gp_mn_zz
    clear G_mn_zz_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{3,3,1}(:,:,:) = fftn(fN_ch_all{3,3,1}(:,:,:));
        disp(['Time for FFT of Gzz_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gzz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_zz_e = ten_mat_prod(fN_charge_Toeplitz{34},{fN_charge_Toeplitz{35},fN_charge_Toeplitz{36},fN_charge_Toeplitz{37}});
        [fN_ch_all{3,3,2}] = Toeplitz_JVIE('G_zz_Efield',L+1,M+1,N+1,dx^2*G_mn_zz_e(1:L,1:M,1:N+1)); % Gp_mn_zz
        clear G_mn_zz_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{3,3,2}(:,:,:) = fftn(fN_ch_all{3,3,2}(:,:,:));
            disp(['Time for FFT of Gzz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_zz_factor_matrix1,d_zz_factor_matrix2,d_zz_factor_matrix3,d_zz_core_tensor,d_zz_mem] = Tucker(fN_ch_all{3,3,2},tolerance);
                fN_ch_all{3,3,2}=[];
                fN_charge{46}=d_zz_core_tensor;% dielectric_zz
                fN_charge{47}=d_zz_factor_matrix1;
                fN_charge{48}=d_zz_factor_matrix2;
                fN_charge{49}=d_zz_factor_matrix3;
                clear d_zz_factor_matrix1 d_zz_factor_matrix2 d_zz_factor_matrix3 d_zz_core_tensor  ;
                disp(['Time for compression of Gzz_d ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{3,3,1}] = Toeplitz_JVIE('G_zz_charge',L+1,M+1,N+1,dx^3.*fN_charge_Toeplitz{3,3,1}(1:L,1:M,1:N+1)); % Gp_mn_zz
    if (fl_no_fft == 0)
        tic
        fN_ch_all{3,3,1}(:,:,:) = fftn(fN_ch_all{3,3,1}(:,:,:));
        disp(['Time for FFT of Gzz_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gzz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        [fN_ch_all{3,3,2}] = Toeplitz_JVIE('G_zz_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{3,3,2}(1:L,1:M,1:N+1)); % Gp_mn_zz
        if (fl_no_fft == 0)
            tic
            fN_ch_all{3,3,2}(:,:,:) = fftn(fN_ch_all{3,3,2}(:,:,:));
            disp(['Time for FFT of Gzz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_zz_factor_matrix1,d_zz_factor_matrix2,d_zz_factor_matrix3,d_zz_core_tensor,d_zz_mem] = Tucker(fN_ch_all{3,3,2},tolerance);
                fN_ch_all{3,3,2}=[];
                fN_charge{46}=d_zz_core_tensor;% dielectric_zz
                fN_charge{47}=d_zz_factor_matrix1;
                fN_charge{48}=d_zz_factor_matrix2;
                fN_charge{49}=d_zz_factor_matrix3;
                clear d_zz_factor_matrix1 d_zz_factor_matrix2 d_zz_factor_matrix3 d_zz_core_tensor  ;
                disp(['Time for compression of Gzz_d ::: ',num2str(toc)])
            end
        end
    end
end

disp(['Total time for computing circulant tensor of Gzz ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gxy interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_xy_c = ten_mat_prod(fN_charge_Toeplitz{14},{fN_charge_Toeplitz{15},fN_charge_Toeplitz{16},fN_charge_Toeplitz{17}});
    [fN_ch_all{1,2,1}] = Toeplitz_JVIE('G_xy_charge',L+1,M+1,N+1,dx^3*G_mn_xy_c(1:L+2,1:M+2,1:N)); % Gp_mn_xy
    clear G_mn_xy_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,2,1}(:,:,:) = fftn(fN_ch_all{1,2,1}(:,:,:));
        disp(['Time for FFT of Gxy_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gxy_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_xy_e = ten_mat_prod(fN_charge_Toeplitz{38},{fN_charge_Toeplitz{39},fN_charge_Toeplitz{40},fN_charge_Toeplitz{41}});
        [fN_ch_all{1,2,2}] = Toeplitz_JVIE('G_xy_Efield',L+1,M+1,N+1,dx^2*G_mn_xy_e(1:L+2,1:M+2,1:N)); % Gp_mn_xy
        clear G_mn_xy_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{1,2,2}(:,:,:) = fftn(fN_ch_all{1,2,2}(:,:,:));
            disp(['Time for FFT of Gxy_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xy_factor_matrix1,d_xy_factor_matrix2,d_xy_factor_matrix3,d_xy_core_tensor,d_xy_mem] = Tucker(fN_ch_all{1,2,2},tolerance);
                fN_ch_all{1,2,2}=[];
                fN_charge{50}=d_xy_core_tensor;% dielectric_xy
                fN_charge{51}=d_xy_factor_matrix1;
                fN_charge{52}=d_xy_factor_matrix2;
                fN_charge{53}=d_xy_factor_matrix3;
                clear d_xy_factor_matrix1 d_xy_factor_matrix2 d_xy_factor_matrix3 d_xy_core_tensor  ;
                disp(['Time for compression of Gxy_d ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{1,2,1}] = Toeplitz_JVIE('G_xy_charge',L+1,M+1,N+1,dx^3.*fN_charge_Toeplitz{1,2,1}(1:L+2,1:M+2,1:N)); % Gp_mn_xy
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,2,1}(:,:,:) = fftn(fN_ch_all{1,2,1}(:,:,:));
        disp(['Time for FFT of Gxy_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gxy_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        [fN_ch_all{1,2,2}] = Toeplitz_JVIE('G_xy_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{1,2,2}(1:L+2,1:M+2,1:N)); % Gp_mn_xy
        if (fl_no_fft == 0)
            tic
            fN_ch_all{1,2,2}(:,:,:) = fftn(fN_ch_all{1,2,2}(:,:,:));
            disp(['Time for FFT of Gxy_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xy_factor_matrix1,d_xy_factor_matrix2,d_xy_factor_matrix3,d_xy_core_tensor,d_xy_mem] = Tucker(fN_ch_all{1,2,2},tolerance);
                fN_ch_all{1,2,2}=[];
                fN_charge{50}=d_xy_core_tensor;% dielectric_xy
                fN_charge{51}=d_xy_factor_matrix1;
                fN_charge{52}=d_xy_factor_matrix2;
                fN_charge{53}=d_xy_factor_matrix3;
                clear d_xy_factor_matrix1 d_xy_factor_matrix2 d_xy_factor_matrix3 d_xy_core_tensor  ;
                disp(['Time for compression of Gxy_d ::: ',num2str(toc)])
            end
        end
    end
end

disp(['TOtal time for computing circulant tensor of Gxy ::: ',num2str(toc)])
%%%---------------------------------------------------
% Gxz interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_xz_c = ten_mat_prod(fN_charge_Toeplitz{18},{fN_charge_Toeplitz{19},fN_charge_Toeplitz{20},fN_charge_Toeplitz{21}});
    [fN_ch_all{1,3,1}] = Toeplitz_JVIE('G_xz_charge',L+1,M+1,N+1,dx^3*G_mn_xz_c(1:L+2,1:M,1:N+2)); % Gp_mn_xz
    clear G_mn_xz_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,3,1}(:,:,:) = fftn(fN_ch_all{1,3,1}(:,:,:));
        disp(['Time for FFT of Gxz_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gxz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_xz_e = ten_mat_prod(fN_charge_Toeplitz{42},{fN_charge_Toeplitz{43},fN_charge_Toeplitz{44},fN_charge_Toeplitz{45}});
        [fN_ch_all{1,3,2}] = Toeplitz_JVIE('G_xz_Efield',L+1,M+1,N+1,dx^2*G_mn_xz_e(1:L+2,1:M,1:N+2)); % Gp_mn_xz
        clear G_mn_xz_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{1,3,2}(:,:,:) = fftn(fN_ch_all{1,3,2}(:,:,:));
            disp(['Time for FFT of Gxz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xz_factor_matrix1,d_xz_factor_matrix2,d_xz_factor_matrix3,d_xz_core_tensor,d_xz_mem] = Tucker(fN_ch_all{1,3,2},tolerance);
                fN_ch_all{1,3,2}=[];
                fN_charge{54}=d_xz_core_tensor;% dielectric_xz
                fN_charge{55}=d_xz_factor_matrix1;
                fN_charge{56}=d_xz_factor_matrix2;
                fN_charge{57}=d_xz_factor_matrix3;
                clear d_xz_factor_matrix1 d_xz_factor_matrix2 d_xz_factor_matrix3 d_xz_core_tensor  ;
                disp(['Time for compression of Gxz ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{1,3,1}] = Toeplitz_JVIE('G_xz_charge',L+1,M+1,N+1,dx^3.*fN_charge_Toeplitz{1,3,1}(1:L+2,1:M,1:N+2)); % Gp_mn_xz
    if (fl_no_fft == 0)
        tic
        fN_ch_all{1,3,1}(:,:,:) = fftn(fN_ch_all{1,3,1}(:,:,:));
        disp(['Time for FFT of Gxz_c ::: ',num2str(toc)])
        if fl_Tucker_decomp==1
            tic
            [LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
            [c_xz_factor_matrix1,c_xz_factor_matrix2,c_xz_factor_matrix3,c_xz_core_tensor,c_xz_mem] = Tucker(fN_ch_all{1,3,1},tolerance);
            [c_zx_factor_matrix1,c_zx_factor_matrix2,c_zx_factor_matrix3,c_zx_core_tensor,c_zx_mem] = Tucker(cong(fN_ch_all{1,3,1}),tolerance);
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
            disp(['Time for compression of Gxz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
         [fN_ch_all{1,3,2}] = Toeplitz_JVIE('G_xz_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{1,3,2}(1:L+2,1:M,1:N+2)); % Gp_mn_xz
         if (fl_no_fft == 0)
            tic
            fN_ch_all{1,3,2}(:,:,:) = fftn(fN_ch_all{1,3,2}(:,:,:));
            disp(['Time for FFT of Gxz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_xz_factor_matrix1,d_xz_factor_matrix2,d_xz_factor_matrix3,d_xz_core_tensor,d_xz_mem] = Tucker(fN_ch_all{1,3,2},tolerance);
                fN_ch_all{1,3,2}=[];
                fN_charge{54}=d_xz_core_tensor;% dielectric_xz
                fN_charge{55}=d_xz_factor_matrix1;
                fN_charge{56}=d_xz_factor_matrix2;
                fN_charge{57}=d_xz_factor_matrix3;
                clear d_xz_factor_matrix1 d_xz_factor_matrix2 d_xz_factor_matrix3 d_xz_core_tensor  ;
                disp(['Time for compression of Gxz ::: ',num2str(toc)])
            end
        end
    end
end
disp(['Total time for computing circulant tensor of Gxz ::: ',num2str(toc)])

%%%---------------------------------------------------
% Gyz interactions
tic
if size_stored_tensor(1) > 3 %the prestored tensor is compressed
    G_mn_yz_c = ten_mat_prod(fN_charge_Toeplitz{22},{fN_charge_Toeplitz{23},fN_charge_Toeplitz{24},fN_charge_Toeplitz{25}});
    [fN_ch_all{2,3,1}] = Toeplitz_JVIE('G_yz_charge',L+1,M+1,N+1,dx^3*G_mn_yz_c(1:L,1:M+2,1:N+2)); % Gp_mn_yz
    clear G_mn_yz_c
    if (fl_no_fft == 0)
        tic
        fN_ch_all{2,3,1}(:,:,:) = fftn(fN_ch_all{2,3,1}(:,:,:));
        disp(['Time for FFT of Gyz_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gyz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
        G_mn_yz_e = ten_mat_prod(fN_charge_Toeplitz{46},{fN_charge_Toeplitz{47},fN_charge_Toeplitz{48},fN_charge_Toeplitz{49}});
        [fN_ch_all{2,3,2}] = Toeplitz_JVIE('G_yz_Efield',L+1,M+1,N+1,dx^2*G_mn_yz_e(1:L,1:M+2,1:N+2)); % Gp_mn_yz
        clear G_mn_yz_e
        if (fl_no_fft == 0)
            tic
            fN_ch_all{2,3,2}(:,:,:) = fftn(fN_ch_all{2,3,2}(:,:,:));
            disp(['Time for FFT of Gyz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_yz_factor_matrix1,d_yz_factor_matrix2,d_yz_factor_matrix3,d_yz_core_tensor,d_yz_mem] = Tucker(fN_ch_all{2,3,2},tolerance);
                fN_ch_all{2,3,2}=[];
                fN_charge{58}=d_yz_core_tensor;% dielectric_yz
                fN_charge{59}=d_yz_factor_matrix1;
                fN_charge{60}=d_yz_factor_matrix2;
                fN_charge{61}=d_yz_factor_matrix3;
                clear d_yz_factor_matrix1 d_yz_factor_matrix2 d_yz_factor_matrix3 d_yz_core_tensor  ;
                disp(['Time for compression of Gyz_d ::: ',num2str(toc)])
            end
        end
    end
else                         %the prestored tensor is not compressed
    [fN_ch_all{2,3,1}] = Toeplitz_JVIE('G_yz_charge',L+1,M+1,N+1,dx^3.*fN_charge_Toeplitz{2,3,1}(1:L,1:M+2,1:N+2)); % Gp_mn_yz
    if (fl_no_fft == 0)
        tic
        fN_ch_all{2,3,1}(:,:,:) = fftn(fN_ch_all{2,3,1}(:,:,:));
        disp(['Time for FFT of Gyz_c ::: ',num2str(toc)])
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
            disp(['Time for compression of Gyz_c ::: ',num2str(toc)])
        end
    end
    if (fl_diel == 1)
         [fN_ch_all{2,3,2}] = Toeplitz_JVIE('G_yz_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{2,3,2}(1:L,1:M+2,1:N+2)); % Gp_mn_yz
         if (fl_no_fft == 0)
            tic
            fN_ch_all{2,3,2}(:,:,:) = fftn(fN_ch_all{2,3,2}(:,:,:));
            disp(['Time for FFT of Gyz_d ::: ',num2str(toc)])
            if fl_Tucker_decomp==1
                tic
                [d_yz_factor_matrix1,d_yz_factor_matrix2,d_yz_factor_matrix3,d_yz_core_tensor,d_yz_mem] = Tucker(fN_ch_all{2,3,2},tolerance);
                fN_ch_all{2,3,2}=[];
                fN_charge{58}=d_yz_core_tensor;% dielectric_yz
                fN_charge{59}=d_yz_factor_matrix1;
                fN_charge{60}=d_yz_factor_matrix2;
                fN_charge{61}=d_yz_factor_matrix3;
                clear d_yz_factor_matrix1 d_yz_factor_matrix2 d_yz_factor_matrix3 d_yz_core_tensor  ;
                disp(['Time for compression of Gyz_d ::: ',num2str(toc)])
            end
        end
    end
end

disp(['Time for computing circulant tensor of Gyz ::: ',num2str(toc)])

%%%---------------------------------------------------
% new ones (involving off-diagonal blocks) here just for dielectric !!!
if (fl_diel == 1)
    % Gyx interactions
    tic
    if size_stored_tensor(1) > 3 %the prestored tensor is compressed
        G_mn_yx_e = ten_mat_prod(fN_charge_Toeplitz{50},{fN_charge_Toeplitz{51},fN_charge_Toeplitz{52},fN_charge_Toeplitz{53}});
        [fN_ch_all{2,1,2}] = Toeplitz_JVIE('G_yx_Efield',L+1,M+1,N+1,dx^2*G_mn_yx_e(1:L+2,1:M+2,1:N)); % Gp_mn_yx
        clear G_mn_yx_e
    else                         %the prestored tensor is not compressed
        [fN_ch_all{2,1,2}] = Toeplitz_JVIE('G_yx_Efield',L+1,M+1,N+1, dx^2*fN_charge_Toeplitz{2,1,2}(1:L+2,1:M+2,1:N)); % Gp_mn_yx
    end
    
    disp(['Time for computing circulant tensor of Gyx ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        fN_ch_all{2,1,2}(:,:,:) = fftn(fN_ch_all{2,1,2}(:,:,:));
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
    end
    %%%---------------------------------------------------
    % Gzx interactions
    tic
    if size_stored_tensor(1) > 3 %the prestored tensor is compressed
        G_mn_zx_e = ten_mat_prod(fN_charge_Toeplitz{54},{fN_charge_Toeplitz{55},fN_charge_Toeplitz{56},fN_charge_Toeplitz{57}});
        [fN_ch_all{3,1,2}] = Toeplitz_JVIE('G_zx_Efield',L+1,M+1,N+1,dx^2*G_mn_zx_e(1:L+2,1:M,1:N+2)); % Gp_mn_zx
        clear G_mn_zx_e
    else                         %the prestored tensor is not compressed
        [fN_ch_all{3,1,2}] = Toeplitz_JVIE('G_zx_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{3,1,2}(1:L+2,1:M,1:N+2)); % Gp_mn_zx
    end
    
    disp(['Time for computing circulant tensor of Gzx ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & cross_Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        fN_ch_all{3,1,2}(:,:,:) = fftn(fN_ch_all{3,1,2}(:,:,:));
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
    end
    %%%---------------------------------------------------
    % Gzy interactions
    tic
    if size_stored_tensor(1) > 3 %the prestored tensor is compressed
        G_mn_zy_e = ten_mat_prod(fN_charge_Toeplitz{58},{fN_charge_Toeplitz{59},fN_charge_Toeplitz{60},fN_charge_Toeplitz{61}});
        [fN_ch_all{3,2,2}] = Toeplitz_JVIE('G_zy_Efield',L+1,M+1,N+1,dx^2*G_mn_zy_e(1:L,1:M+2,1:N+2)); % Gp_mn_zx
        clear G_mn_zy_e
    else                         %the prestored tensor is not compressed
        [fN_ch_all{3,2,2}] = Toeplitz_JVIE('G_zy_Efield',L+1,M+1,N+1,dx^2*fN_charge_Toeplitz{3,2,2}(1:L,1:M+2,1:N+2)); % Gp_mn_zy
    end
    
    disp(['Time for computing circulant tensor of Gzy ::: ',num2str(toc)])
    %%%----------------------------------------------
    %               FFT & cross_Tucker
    %%%----------------------------------------------
    if (fl_no_fft == 0)
        tic
        fN_ch_all{3,2,2}(:,:,:) = fftn(fN_ch_all{3,2,2}(:,:,:));
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
disp('Done... Computing circulant tensors    ')
disp('-----------------------------------------------------')