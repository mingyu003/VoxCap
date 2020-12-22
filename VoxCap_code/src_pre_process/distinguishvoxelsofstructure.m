function [idx,grid_intcon] = distinguishvoxelsofstructure(r,Res,Cnt,Dims,Orients,Eps_inout,fl_plot_vox_structure)
%%    Defines constitutive parameters of voxels and 
%     distinguishes which voxel belongs to which dielectric or conductor

%% INPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%   Cnt         Cartesian coordinates of the centers of
%               conds/diels
%               row: cond/diel id, column: dimensions
%   Dims        Dimensions of conds/diels(LxWxH) (along x, y, and z)
%               row: cond/diel id, column: length,width,height
%   Orients     Orientations of cond/diel, string
%   Eps_inout   Relative permittivity of inner and outer media of
%               dielectric and conductor - Assumption: epsr of inner medium
%               of conductor is assumed to be 0. 
%               row: cond/diel id, column: inner rel. perm., outer rel. perm.
%
%% OUTPUT
%   grid_intcon  4D (LxMxNx3) array with domain voxelized structure
%                coordinates; the coordinates(0,0,0) is assigned to air
%                voxels. Note that this is dangerous!
%   idx          Vector storing the ids of non-zero elements in grid_intcon
%   idx_c_or_d   Vector storing 1 for conductor and 0 for dielectric 
%                for the elements in idx
%   idx_cd_no    Vector storing the ids of conductor or dielectric structures 
%                for the elements in idx

%
% -------------------------------------------------------------------------
%
%   Abdulkadir C. Yucel -- acyucel@ntu.edu.sg
%   Computational Electromagnetics Group, NTU
%
% _________________________________________________________________________

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);

num_indv_elem=size(Cnt,1);

% define bounds for each interconnect
x_bnd=zeros(num_indv_elem,2);y_bnd=zeros(num_indv_elem,2);z_bnd=zeros(num_indv_elem,2);
for kk=1:num_indv_elem
    temp_cen=Cnt(kk,1:3);
    temp_dim=Dims(kk,1:3);
    
    if (Orients(kk) == 'x')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(1)*0.5 temp_cen(1)+temp_dim(1)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5];
    elseif (Orients(kk) == 'y')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(2)*0.5 temp_cen(1)+temp_dim(2)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(1)*0.5 temp_cen(2)+temp_dim(1)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5];
    elseif (Orients(kk) == 'z')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(3)*0.5 temp_cen(1)+temp_dim(3)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(1)*0.5 temp_cen(3)+temp_dim(1)*0.5];
    else
        error('Orients should have only x, y, or z strings')
    end
    
end

boolean_tens=zeros(L,M,N,num_indv_elem);

tola=1e-12;
for ll=1:L
    for mm=1:M
        for nn=1:N
            for kk=1:num_indv_elem
                if ( r(ll,mm,nn,1) > x_bnd(kk,1)-tola && r(ll,mm,nn,1) < x_bnd(kk,2)+tola &&...
                        r(ll,mm,nn,2) > y_bnd(kk,1)-tola && r(ll,mm,nn,2) < y_bnd(kk,2)+tola && ...
                        r(ll,mm,nn,3) > z_bnd(kk,1)-tola && r(ll,mm,nn,3) < z_bnd(kk,2)+tola)
                    
                    boolean_tens(ll,mm,nn,kk)=1;

                end
                
            end
        end
    end
end

idx=cell(num_indv_elem,1);
for kk=1:num_indv_elem
    idx{kk} = find(boolean_tens(:,:,:,kk)==1);  % get indices of elements
end

% st_is_c_or_d=zeros(num_indv_elem,1); % for conductor ->1, for dielectric->0
% for kk=1:num_indv_elem
%     if (abs(Eps_inout(kk,1))<0)
%       st_is_c_or_d(kk)=1;  
%     end
% end


% -------------------------------------------------------------------------
% Assign the grid of conductors and visualize it 
% -------------------------------------------------------------------------

grid_intcon=cell(num_indv_elem,1);
for kk=1:num_indv_elem
    [L,M,N,~] = size(r);
    grid_intcon{kk} = zeros(L,M,N,3);
    grid_intcon{kk}(idx{kk}) = r(idx{kk});
    grid_intcon{kk}(L*M*N+idx{kk}) = r(L*M*N+idx{kk});
    grid_intcon{kk}(2*L*M*N+idx{kk}) = r(2*L*M*N+idx{kk});
    
    % disp('Attention ::: Pad NaNs to grid_intcon for air voxels instead zeros! ')
    % disp('This will avoid possible bug that appears when you have a voxel centered at (0,0,0)')
    
    % plot geometry
    %fl_plot_vox_structure=1;
    if (fl_plot_vox_structure == 1)
        tic
        plot_boxes_of_grid(grid_intcon{kk},Res);
        title('Voxelized Individual Structure')
        disp(['Time for visualizing the seperate body in voxelized structure ::: ',num2str(toc)])
    end
end
