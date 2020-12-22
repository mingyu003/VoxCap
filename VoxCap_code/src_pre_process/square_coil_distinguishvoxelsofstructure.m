function [idx,grid_intcon,boolean_tens] = square_coil_distinguishvoxelsofstructure(dx,r,Cnt,Dims,Orients,num_cond_x,num_cond_z)
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
tic;
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
disp(['Time for xyz_bnd ::: ',num2str(toc)])

tic;
idx=cell(num_indv_elem,1);
% tola=1e-12;

for kk=1:num_indv_elem
%     counter=0;
%     idx{kk}=zeros(L*M*N,1);
    boolean_tens=zeros(L*M*N,1);
    %      boolean_tens=boolean_loop(x_bnd(kk,:),y_bnd(kk,:),z_bnd(kk,:),r);
    for ll=round((x_bnd(kk,1))/dx+1):round((x_bnd(kk,2))/dx)
        for mm=round((y_bnd(kk,1))/dx+1):round((y_bnd(kk,2))/dx)
            for nn=round((z_bnd(kk,1))/dx+1):round((z_bnd(kk,2))/dx)
                
                if ( r(ll,mm,nn,1) > x_bnd(kk,1)-tola && r(ll,mm,nn,1) < x_bnd(kk,2)+tola &&...
                        r(ll,mm,nn,2) > y_bnd(kk,1)-tola && r(ll,mm,nn,2) < y_bnd(kk,2)+tola && ...
                        r(ll,mm,nn,3) > z_bnd(kk,1)-tola && r(ll,mm,nn,3) < z_bnd(kk,2)+tola) 
                    boolean_tens(ll+(mm-1)*L+(nn-1)*L*M)=1;
%                     idx{kk}=[idx{kk};ll+(mm-1)*L+(nn-1)*L*M];
%                     counter=counter+1;
                end
            end
        end
    end  
%     idx{kk}=idx{kk}(1:counter);
    idx{kk} = find(boolean_tens(:)==1);  % get indices of elements 
end
clear boolean_tens x_bnd y_bnd z_bnd; 
disp(['Time for boolean loop ::: ',num2str(toc)])
tic
idx_new=cell(num_cond_z,1);
ids=cell(num_cond_z,1);

for ii=1:num_cond_z
    for jj=1:2*num_cond_x
        ids{ii}=find(Cnt(:,3)==Cnt(ii,3));
        idx_new{ii}=[idx_new{ii};idx{ids{ii}(jj)}];
    end
end
idx=idx_new;
clear idx_new;
disp(['Time for new idx ::: ',num2str(toc)])
% -------------------------------------------------------------------------
% Assign the grid of conductors and visualize it 
% -------------------------------------------------------------------------
tic;
grid_intcon=cell(num_cond_z,1);
for kk=1:num_cond_z
%     [L,M,N,~] = size(r);
    grid_intcon{kk} = zeros(L,M,N,3);
    grid_intcon{kk}(idx{kk}) = r(idx{kk});
    grid_intcon{kk}(L*M*N+idx{kk}) = r(L*M*N+idx{kk});
    grid_intcon{kk}(2*L*M*N+idx{kk}) = r(2*L*M*N+idx{kk});
end
disp(['Time for grid ::: ',num2str(toc)])
