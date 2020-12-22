function [Gp_mn] = precond_Toeplitz(fl_block,Nx,Ny,Nz,G_mn)


switch fl_block
    
    case 'G_xx_yy_zz_charge'
        
        % Nx->L, Ny->M, Nz->N
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn;
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)         = G_mn(Nx:-1:2,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)         = G_mn(1:Nx,Ny:-1:2,1:Nz);
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)         = G_mn(1:Nx,1:Ny,Nz:-1:2);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,1:Nz)     = G_mn(Nx:-1:2,Ny:-1:2,1:Nz);
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,Nz+2:2*Nz)     = G_mn(Nx:-1:2,1:Ny,Nz:-1:2);
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,Nz+2:2*Nz)     = G_mn(1:Nx,Ny:-1:2,Nz:-1:2);
        % Cube 'LMN'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,Nz+2:2*Nz) = G_mn(Nx:-1:2,Ny:-1:2,Nz:-1:2);
        
        
    case 'G_xx_Efield' % only x-variation or blocks related to L will have negative sign
        
        % Nx->L, Ny->M, Nz->N
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn;
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)         = -G_mn(Nx:-1:2,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)         = G_mn(1:Nx,Ny:-1:2,1:Nz);
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)         = G_mn(1:Nx,1:Ny,Nz:-1:2);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,1:Nz)     = -G_mn(Nx:-1:2,Ny:-1:2,1:Nz);
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,Nz+2:2*Nz)     = -G_mn(Nx:-1:2,1:Ny,Nz:-1:2);
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,Nz+2:2*Nz)     = G_mn(1:Nx,Ny:-1:2,Nz:-1:2);
        % Cube 'LMN'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,Nz+2:2*Nz) = -G_mn(Nx:-1:2,Ny:-1:2,Nz:-1:2);
        
    case 'G_yy_Efield' % only y-variation or blocks related to M will have negative sign
        
        % Nx->L, Ny->M, Nz->N
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn;
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)         = G_mn(Nx:-1:2,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)         = -G_mn(1:Nx,Ny:-1:2,1:Nz);
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)         = G_mn(1:Nx,1:Ny,Nz:-1:2);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,1:Nz)     = -G_mn(Nx:-1:2,Ny:-1:2,1:Nz);
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,Nz+2:2*Nz)     = G_mn(Nx:-1:2,1:Ny,Nz:-1:2);
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,Nz+2:2*Nz)     = -G_mn(1:Nx,Ny:-1:2,Nz:-1:2);
        % Cube 'LMN'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,Nz+2:2*Nz) = -G_mn(Nx:-1:2,Ny:-1:2,Nz:-1:2);
        
    case 'G_zz_Efield' % only z-variation or blocks related to N will have negative sign
        
        % Nx->L, Ny->M, Nz->N
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn;
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)         = G_mn(Nx:-1:2,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)         = G_mn(1:Nx,Ny:-1:2,1:Nz);
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)         = -G_mn(1:Nx,1:Ny,Nz:-1:2);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,1:Nz)     = G_mn(Nx:-1:2,Ny:-1:2,1:Nz);
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,Nz+2:2*Nz)     = -G_mn(Nx:-1:2,1:Ny,Nz:-1:2);
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,Nz+2:2*Nz)     = -G_mn(1:Nx,Ny:-1:2,Nz:-1:2);
        % Cube 'LMN'
        Gp_mn(Nx+2:2*Nx,Ny+2:2*Ny,Nz+2:2*Nz) = -G_mn(Nx:-1:2,Ny:-1:2,Nz:-1:2);
        
        
    case 'G_xy_charge'
        
        % Nx->L+1, Ny->M+1, Nz->N
        % The dimensions of G_xy Toeplitz are (L+2)x(M+2)x(N)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,1:Nz) =  G_mn(1:Nx,1:1:Ny-1,1:Nz);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,2*Ny:-1:Ny+2,1:Nz)  =  G_mn(Nx+1:-1:3,1:Ny-1,1:Nz);
        % Cube 'N', 'LM', 'LN', 'LMN' all together
        Gp_mn(:,:,2*Nz:-1:Nz+2)  =  Gp_mn(:,:,2:Nz);
        
        
    case 'G_xy_Efield' % only x-variation or blocks related to L will have negative sign [as derivative is on observer x-directed panel]
        
        % Nx->L+1, Ny->M+1, Nz->N
        % The dimensions of G_xy Toeplitz are (L+2)x(M+2)x(N)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        % Cube 'M'
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,1:Nz) =  G_mn(1:Nx,1:1:Ny-1,1:Nz);
        % Cube 'LM'
        Gp_mn(Nx+2:2*Nx,2*Ny:-1:Ny+2,1:Nz)  =  -G_mn(Nx+1:-1:3,1:Ny-1,1:Nz);
        % Cube 'N', 'LM', 'LN', 'LMN' all together
        Gp_mn(:,:,2*Nz:-1:Nz+2)  =  Gp_mn(:,:,2:Nz);
        
    case 'G_yx_charge' % only y-variation or blocks related to M will have negative sign [as derivative is on observer y-directed panel]
        
        % Nx->L+1, Ny->M+1, Nz->N
        % The dimensions of G_xy Toeplitz are (L+2)x(M+2)x(N)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        Gp_mn(2*Nx:-1:Nx+2,1:Ny,1:Nz) =  G_mn(1:Nx-1,1:Ny,1:Nz);
        
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  G_mn(1:Nx,Ny+1:-1:3,1:Nz);
        
        % Cube 'LM'
        Gp_mn(2*Nx:-1:Nx+2,Ny+2:2*Ny,1:Nz)  =  G_mn(1:Nx-1,Ny+1:-1:3,1:Nz);
        
        % Cube 'N', 'LM', 'LN', 'LMN' all together
        Gp_mn(:,:,2*Nz:-1:Nz+2)  =  Gp_mn(:,:,2:Nz);
        
    case 'G_yx_Efield' % only y-variation or blocks related to M will have negative sign [as derivative is on observer y-directed panel]
        
        % Nx->L+1, Ny->M+1, Nz->N
        % The dimensions of G_xy Toeplitz are (L+2)x(M+2)x(N)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
         Gp_mn(2*Nx:-1:Nx+2,1:Ny,1:Nz) =  G_mn(1:Nx-1,1:Ny,1:Nz);
 
         % Cube 'M'
         Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz);
         
         % Cube 'LM'
         Gp_mn(2*Nx:-1:Nx+2,Ny+2:2*Ny,1:Nz)  =  -G_mn(1:Nx-1,Ny+1:-1:3,1:Nz);

         % Cube 'N', 'LM', 'LN', 'LMN' all together
         Gp_mn(:,:,2*Nz:-1:Nz+2)  =  Gp_mn(:,:,2:Nz);        

       
    case 'G_xz_charge'
        
        % Nx->L+1; Ny->M; Nz->N+1;
        % The dimensions of G_xz Toeplitz are (L+2)xMx(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,2*Nz:-1:Nz+2)  =  G_mn(Nx+1:-1:3,1:Ny,1:Nz-1);
        
        % Cube 'M', 'LM', 'MN', 'LMN' all together
        Gp_mn(:,2*Ny:-1:Ny+2,:)  =  Gp_mn(:,2:Ny,:);
       
    case 'G_zx_charge' % only z-variation or blocks related to N will have negative sign [as derivative is on observer z-directed panel]
        
        % Nx->L+1; Ny->M; Nz->N+1;
        % The dimensions of G_xz Toeplitz are (L+2)xMx(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        %Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        Gp_mn(2*Nx:-1:Nx+2,1:Ny,1:Nz) =  G_mn(1:1:Nx-1,1:1:Ny,1:1:Nz);
        
        % Cube 'N'
        %Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)  =  G_mn(1:Nx,1:Ny,Nz+1:-1:3);
        
        % Cube 'LN'
        %Gp_mn(Nx+2:2*Nx,1:Ny,2*Nz:-1:Nz+2)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz-1);
        Gp_mn(2*Nx:-1:Nx+2,1:Ny,Nz+2:2*Nz)  =  G_mn(1:Nx-1,1:Ny,Nz+1:-1:3);
        
        % Cube 'M', 'LM', 'MN', 'LMN' all together
        Gp_mn(:,2*Ny:-1:Ny+2,:)  =  Gp_mn(:,2:Ny,:);
        
    case 'G_xz_Efield'
        
        % Nx->L+1; Ny->M; Nz->N+1;
        % The dimensions of G_xz Toeplitz are (L+2)xMx(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        
        % Cube 'LN'
        Gp_mn(Nx+2:2*Nx,1:Ny,2*Nz:-1:Nz+2)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz-1);
        
        % Cube 'M', 'LM', 'MN', 'LMN' all together
        Gp_mn(:,2*Ny:-1:Ny+2,:)  =  Gp_mn(:,2:Ny,:);
        
        
    case 'G_zx_Efield' % only z-variation or blocks related to N will have negative sign [as derivative is on observer z-directed panel]
        
        % Nx->L+1; Ny->M; Nz->N+1;
        % The dimensions of G_xz Toeplitz are (L+2)xMx(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'L'
        %Gp_mn(Nx+2:2*Nx,1:Ny,1:Nz)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz);
        Gp_mn(2*Nx:-1:Nx+2,1:Ny,1:Nz) =  G_mn(1:1:Nx-1,1:1:Ny,1:1:Nz);
        
        % Cube 'N'
        %Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)  =  -G_mn(1:Nx,1:Ny,Nz+1:-1:3);
        
        % Cube 'LN'
        %Gp_mn(Nx+2:2*Nx,1:Ny,2*Nz:-1:Nz+2)  =  -G_mn(Nx+1:-1:3,1:Ny,1:Nz-1);
        Gp_mn(2*Nx:-1:Nx+2,1:Ny,Nz+2:2*Nz)  =  -G_mn(1:Nx-1,1:Ny,Nz+1:-1:3);
        
        % Cube 'M', 'LM', 'MN', 'LMN' all together
        Gp_mn(:,2*Ny:-1:Ny+2,:)  =  Gp_mn(:,2:Ny,:);
        
    case 'G_yz_charge'
        
        % Nx->L; Ny->M+1; Nz->N+1;
        % The dimensions of G_yz Toeplitz are Lx(M+2)x(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  G_mn(1:Nx,Ny+1:-1:3,1:Nz);
        
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,2*Nz:-1:Nz+2)  =  G_mn(1:Nx,Ny+1:-1:3,1:Nz-1);
        
        % Cube 'L', 'LM', 'LN', 'LMN' all together
        Gp_mn(2*Nx:-1:Nx+2,:,:)  =  Gp_mn(2:Nx,:,:);
    
    case 'G_zy_charge'
        
        % Nx->L; Ny->M+1; Nz->N+1;
        % The dimensions of G_yz Toeplitz are Lx(M+2)x(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'M'
        %Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz);
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,1:Nz) =  G_mn(1:Nx,1:1:Ny-1,1:1:Nz);
        
        % Cube 'N'
        %Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)  =  G_mn(1:Nx,1:Ny,Nz+1:-1:3);
        
        % Cube 'MN'
        %Gp_mn(1:Nx,Ny+2:2*Ny,2*Nz:-1:Nz+2)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz-1);
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,Nz+2:2*Nz)  =  G_mn(1:Nx,1:Ny-1,Nz+1:-1:3);
        
        % Cube 'L', 'LM', 'LN', 'LMN' all together
        Gp_mn(2*Nx:-1:Nx+2,:,:)  =  Gp_mn(2:Nx,:,:);
        
    case 'G_yz_Efield'
        
        % Nx->L; Ny->M+1; Nz->N+1;
        % The dimensions of G_yz Toeplitz are Lx(M+2)x(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'M'
        Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz);
        
        % Cube 'N'
        Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        
        % Cube 'MN'
        Gp_mn(1:Nx,Ny+2:2*Ny,2*Nz:-1:Nz+2)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz-1);
        
        % Cube 'L', 'LM', 'LN', 'LMN' all together
        Gp_mn(2*Nx:-1:Nx+2,:,:)  =  Gp_mn(2:Nx,:,:);
        
        
    case 'G_zy_Efield'
        
        % Nx->L; Ny->M+1; Nz->N+1;
        % The dimensions of G_yz Toeplitz are Lx(M+2)x(N+2)
        
        Gp_mn = zeros(2*Nx,2*Ny,2*Nz);
        
        Gp_mn(1:Nx,1:Ny,1:Nz) = G_mn(1:Nx,1:Ny,1:Nz);
        
        % Cube 'M'
        %Gp_mn(1:Nx,Ny+2:2*Ny,1:Nz)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz);
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,1:Nz) =  G_mn(1:Nx,1:1:Ny-1,1:1:Nz);
        
        % Cube 'N'
        %Gp_mn(1:Nx,1:Ny,2*Nz:-1:Nz+2) =  G_mn(1:Nx,1:1:Ny,1:1:Nz-1);
        Gp_mn(1:Nx,1:Ny,Nz+2:2*Nz)  =  -G_mn(1:Nx,1:Ny,Nz+1:-1:3);
        
        % Cube 'MN'
        %Gp_mn(1:Nx,Ny+2:2*Ny,2*Nz:-1:Nz+2)  =  -G_mn(1:Nx,Ny+1:-1:3,1:Nz-1);
        Gp_mn(1:Nx,2*Ny:-1:Ny+2,Nz+2:2*Nz)  =  -G_mn(1:Nx,1:Ny-1,Nz+1:-1:3);
        
        % Cube 'L', 'LM', 'LN', 'LMN' all together
        Gp_mn(2*Nx:-1:Nx+2,:,:)  =  Gp_mn(2:Nx,:,:);
        
        
        
    otherwise
        
        error('Invalid Toeplitz block selection for embedding in a circulant!')
        
end



