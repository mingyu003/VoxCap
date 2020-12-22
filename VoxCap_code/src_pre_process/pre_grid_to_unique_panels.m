function [unique_panel_locs,unique_panel_bndry_locs]=pre_grid_to_unique_panels(dx,x_bnd,y_bnd,z_bnd,grid_intcon,idxS,idxS_cond,fl_check_pnls)

%   Inputs:
%   dx          resolution
%   x_bnd       minimum and maximum values of x in bounding box
%   y_bnd       minimum and maximum values of y in bounding box
%   z_bnd       minimum and maximum values of z in bounding box
%   grid_intcon 4D (LxMxNx3) array with domain voxelized grid coordinates
%   idxS        ids of nonzero entries in grid_intcon   
%   idxS_cond   ids of conductors for each panel  - each panel belongs to which cond 

%   Outputs:
%   unique_panel_locs/unique_panel_bndry_locs stores the locations,
%   orientations, and global ids of unique panels on voxelized geometry and
%   on its boundaries, respectively
%   (:,1:3) the locations
%   (:,4) the orientations: 1-> x-aligned, 2->y-aligned, 3-> z-aligned; 
%   (:,5) global id of panels among the panels of computational domain
%   (:,6) conductor number - ID of the conductor to which the panel belongs
%   (This will be used to form connection between the panels on global grid)


% Test inputs for stand alone execution of this subroutine
% clc; close all; clear all; 
% pre_define_the_path_for_folders
% fl_check_pnls=2;

% % test 1 - one domain conductor
% L=4; M=4; N=1; dx = 0.1; 
% grid_tmp = ones(L,M,N); idxS = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp;
% idxS_cond=ones(length(idxS),1);idxS_cond(length(idxS)/2+1:end)=2;
% idxS_cond_dum=[idxS_cond;idxS_cond;idxS_cond];
% x_bnd=[0 L*dx]; y_bnd=[0 M*dx]; z_bnd=[0 N*dx];
% [grid_intcon] = generategridfrombbox(dx,x_bnd,y_bnd,z_bnd,1);

% % test 2 - two conductors - x aligned
% L=4; M=6; N=2; dx = 0.1; %grid_tmp = ones(L,M,N);
% grid_tmp = zeros(L,M,N); grid_tmp(:,1:2,1)=1; grid_tmp(:,M-1:M,N)=1;
% idxS = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp;
% %idxS_cond=ones(length(idxS),1);idxS_cond(length(idxS)/2+1:end)=2;
% idxS_cond=ones(length(idxS),1);
% idxS_cond(1:4)=2;idxS_cond(5:8)=1;idxS_cond(9:12)=1;idxS_cond(13:16)=2;
% idxS_cond_dum=[idxS_cond;idxS_cond;idxS_cond];
% x_bnd=[0 L*dx]; y_bnd=[0 M*dx]; z_bnd=[0 N*dx];
% [rr] = generategridfrombbox(dx,x_bnd,y_bnd,z_bnd,1);
% grid_intcon = zeros(L,M,N,3);grid_intcon(idxS) = rr(idxS);
% grid_intcon(L*M*N+idxS) = rr(L*M*N+idxS);
% grid_intcon(2*L*M*N+idxS) = rr(2*L*M*N+idxS);

% End of test inputs

fl_profile = 0;
idxS_cond_dum=[idxS_cond;idxS_cond;idxS_cond];

% definition
num_nonair_cube=length(idxS);
[L,M,N,~]=size(grid_intcon);

disp('-----------------------------------------------------')
disp('Extracting unique panels of a grid')
disp('This routine uses graph toolbox of Matlab')

% 1) Get the boolean tensor
tic
boolean_tens=zeros(L,M,N);
boolean_tens(idxS)=1;
if(fl_profile == 1); disp(['Time for obtaining boolean matrix ::: ', num2str(toc)]); end

% 2) Number the non-empty voxels (or nodes in graph) along x, y, and z
% directions. Three different number is used to ensure the correct ordering
% of panels along x, y, and z directions in the graphs.

tic
% 2a) ijk to ind tensor and ind to ijk tensor for numbering along x(x_num_vect)
unkids_ijk_to_ind=zeros(L,M,N);
x_num_vect=zeros(num_nonair_cube,3);

dum=0;
for mm=1:N
    for ll=1:M
        for kk=1:L
            if (boolean_tens(kk,ll,mm) ~=0)
                dum=dum+1;
                unkids_ijk_to_ind(kk,ll,mm)=dum;
                x_num_vect(dum,1:3)=[kk ll mm];
            end
        end
    end
end
clear unkids_ijk_to_ind

%2b) ind to ijk tensor for numbering along y(y_num_vect)
y_num_vect=zeros(num_nonair_cube,3);
dum=0;
for kk=1:L
    for mm=1:N
        for ll=1:M
            if (boolean_tens(kk,ll,mm) ~=0)
                dum=dum+1;
                y_num_vect(dum,1:3)=[kk ll mm];
            end
        end
    end
end

%2c) ind to ijk tensor for numbering along z(z_num_vect)
z_num_vect=zeros(num_nonair_cube,3);
dum=0;
for ll=1:M
    for kk=1:L
        for mm=1:N
            if (boolean_tens(kk,ll,mm) ~=0)
                dum=dum+1;
                z_num_vect(dum,1:3)=[kk ll mm];
            end
        end
    end
end

% 2d) ijk tensor storing the numbering along x, y, and z directions
numbering_x_y_z=zeros(L,M,N,3);
for kk=1:size(x_num_vect,1)
    numbering_x_y_z(x_num_vect(kk,1),x_num_vect(kk,2),x_num_vect(kk,3),1)=kk; % x numbering
    numbering_x_y_z(y_num_vect(kk,1),y_num_vect(kk,2),y_num_vect(kk,3),2)=kk; % y numbering
    numbering_x_y_z(z_num_vect(kk,1),z_num_vect(kk,2),z_num_vect(kk,3),3)=kk; % z numbering
end

%2e) mapping from y and z numbering to x numbering

dum=0;
map_from_y2x_numbering=zeros(num_nonair_cube,1);
map_from_z2x_numbering=zeros(num_nonair_cube,1);
for kk=1:L
    for mm=1:N
        for ll=1:M
            if (boolean_tens(kk,ll,mm) ~=0)
                dum=dum+1;
                map_from_y2x_numbering(dum)=numbering_x_y_z(y_num_vect(dum,1),y_num_vect(dum,2),y_num_vect(dum,3),1);
                map_from_z2x_numbering(dum)=numbering_x_y_z(z_num_vect(dum,1),z_num_vect(dum,2),z_num_vect(dum,3),1);
            end
        end
    end
end

if(fl_profile == 1); disp(['Time for getting numbering matrices ::: ', num2str(toc)]);end
clear x_num_vect y_num_vect z_num_vect


% 3) Finding the near neighbors of each non-empty voxel
tic
temp_crit=1; % box_diff
tola=1e-12;
nn_nums_tens=zeros(L,M,N,3); % number of near neighbors of each voxel
nn_x_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> left, 2-> right neighbors
nn_y_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> front, 2-> back neighbors
nn_z_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> bottom, 2-> up neighbors

for mm=1:N%1 % z variation
    for ll=1:M % y variation
        for kk=1:L % x variation
            if (boolean_tens(kk,ll,mm) ~=0) %non-air voxel (source)
                
                lowbnd_indi=kk-temp_crit;
                upbnd_indi=kk+temp_crit;
                
                lowbnd_indj=ll-temp_crit;
                upbnd_indj=ll+temp_crit;
                
                lowbnd_indk=mm-temp_crit;
                upbnd_indk=mm+temp_crit;
                
                dum=0; % counter for near neighbors of each voxel with index (kk,ll,mm)
                % counters for near neighbor along x,y,and z directions
                dum_nn_x_cnt=0; dum_nn_y_cnt=0; dum_nn_z_cnt=0;
                for cc=lowbnd_indk:upbnd_indk % search along z direction
                    
                    if (cc >= 1 && cc <= N)
                        
                        for bb=lowbnd_indj:upbnd_indj % search along y direction
                            
                            if (bb >= 1 && bb <= M)
                                
                                for aa=lowbnd_indi:upbnd_indi % search along x direction
                                    
                                    if (aa >= 1 && aa<= L )
                                        
                                        if (boolean_tens(aa,bb,cc)==1) %non-air voxel (observer)
                                            
                                            diffind=abs(aa-kk)+abs(bb-ll)+abs(cc-mm);
                                            
                                            if( diffind < 2-tola && diffind > 0+tola)% corner and self neighbors
                                                
                                                if (abs(abs(aa-kk)-1) < tola) % x neighbor
                                                    if (sign(aa-kk) < 0) % left neighbor
                                                        nn_x_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,1),1) = numbering_x_y_z(aa,bb,cc,1);
                                                    else % right neighbor
                                                        nn_x_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,1),2) = numbering_x_y_z(aa,bb,cc,1);
                                                    end
                                                    dum_nn_x_cnt=dum_nn_x_cnt+1;
                                                elseif (abs(abs(bb-ll)-1) < tola) % y neighbor
                                                    if (sign(bb-ll) < 0) % front neighbor
                                                        nn_y_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,2),1) = numbering_x_y_z(aa,bb,cc,2);
                                                    else % back neighbor
                                                        nn_y_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,2),2) = numbering_x_y_z(aa,bb,cc,2);
                                                    end
                                                    dum_nn_y_cnt=dum_nn_y_cnt+1;
                                                elseif (abs(abs(cc-mm)-1) < tola) % z neighbor
                                                    if (sign(cc-mm) < 0) % bottom neighbor
                                                        nn_z_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,3),1) = numbering_x_y_z(aa,bb,cc,3);
                                                    else % up neighbor
                                                        nn_z_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,3),2) = numbering_x_y_z(aa,bb,cc,3);
                                                    end
                                                    dum_nn_z_cnt=dum_nn_z_cnt+1;
                                                end
                                                dum=dum+1;
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                nn_nums_tens(kk,ll,mm,1)=dum_nn_x_cnt;
                nn_nums_tens(kk,ll,mm,2)=dum_nn_y_cnt;
                nn_nums_tens(kk,ll,mm,3)=dum_nn_z_cnt;
                
            end
            
        end
        
    end
    
end

clear boolean_tens

if(fl_profile == 1); disp(['Time for determining neighbors of voxels ::: ', num2str(toc)]); end

% 4) Create adjacency matrix for the nodes with numbering along x, y,and
% z directions
tic
% Note: if we create the indices for sparse adjacency matrix and then
% allocate the sparse matrix with sparse command, we see that matrix
% doesn't have the dimension of (num_nonair_cubexnum_nonair_cube)
% (if the last elements are zeros) We should specify the dimensions of
% sparse matrices in sparse command

% adjacency matrix for x directed panels
inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,1)))),3);
dum=0;
for kk=1:size(nn_x_neigh_ind_tens,1)
    if(nn_x_neigh_ind_tens(kk,1) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_x_neigh_ind_tens(kk,1) 1];
    end
    if(nn_x_neigh_ind_tens(kk,2) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_x_neigh_ind_tens(kk,2) 1];
    end
end
clear nn_x_neigh_ind_tens
adj_mat_x=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);

if(isempty(adj_mat_x) == 1)
    adj_mat_x=sparse(num_nonair_cube,num_nonair_cube);
end

% adjacency matrix for y directed panels

inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,2)))),3);
dum=0;
for kk=1:size(nn_y_neigh_ind_tens,1)
    if(nn_y_neigh_ind_tens(kk,1) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_y_neigh_ind_tens(kk,1) 1];
    end
    if(nn_y_neigh_ind_tens(kk,2) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_y_neigh_ind_tens(kk,2) 1];
    end
end
clear nn_y_neigh_ind_tens

adj_mat_y=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);

if(isempty(adj_mat_y) == 1)
    adj_mat_y=sparse(num_nonair_cube,num_nonair_cube);
end

% adjacency matrix for z directed panels

inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,3)))),3);
dum=0;
for kk=1:size(nn_z_neigh_ind_tens,1)
    if(nn_z_neigh_ind_tens(kk,1) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_z_neigh_ind_tens(kk,1) 1];
    end
    if(nn_z_neigh_ind_tens(kk,2) > 0)
        dum=dum+1;
        inds_dum(dum,1:3)=[kk nn_z_neigh_ind_tens(kk,2) 1];
    end
end
clear nn_z_neigh_ind_tens nn_nums_tens

adj_mat_z=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);

if(isempty(adj_mat_z) == 1)
    adj_mat_z=sparse(num_nonair_cube,num_nonair_cube);
end

if(fl_profile == 1); disp(['Time for getting adjacency matrices ::: ', num2str(toc)]); end;

%disp('Are adjacency matrices for x,y,and z symmetric?')
%[issymmetric(adj_mat_x) issymmetric(adj_mat_y) issymmetric(adj_mat_z)]


% 5) Create graphs using adjacency matrices
tic
G_x = graph(adj_mat_x);
G_y = graph(adj_mat_y);
G_z = graph(adj_mat_z);

clear adj_mat_x adj_mat_y adj_mat_z
if(fl_profile == 1); disp(['Time for creating graphs ::: ', num2str(toc)]); end;

% 6) Find the connected elements in each graph
tic
bins_x=conncomp(G_x,'OutputForm','cell');
bins_y=conncomp(G_y,'OutputForm','cell');
bins_z=conncomp(G_z,'OutputForm','cell');
if(fl_profile == 1); disp(['Time for finding conn elements in graphs ::: ', num2str(toc)]);end;

% 7) Find the panels of each voxel (due to x numbering)
tic
% Numbering panels along x direction
ind_panel=0; % counter for all panels along x
voxel2panel_x=zeros(num_nonair_cube,2); % 1 left, 2 right panel
for kk=1:length(bins_x)
    
    for ll=1:length(bins_x{kk})
        
        if (ll == 1) % beginning of connected elements, put a panel
            ind_panel=ind_panel+1;
        end
        
        unk_id=bins_x{kk}(ll);
        
        voxel2panel_x(unk_id,1)=ind_panel; % left panel
        
        ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
        
        voxel2panel_x(unk_id,2)=ind_panel; % right panel
    end
    
end
clear bins_x
num_panel_x=max(max(voxel2panel_x));

% Numbering panels along y direction
ind_panel=0; % counter for all panels along y
voxel2panel_y=zeros(num_nonair_cube,2); % 1 front, 2 back
for kk=1:length(bins_y)
    
    for ll=1:length(bins_y{kk})
        
        if (ll == 1) % beginning of connected elements, put a panel
            ind_panel=ind_panel+1;
        end
        
        unk_id=map_from_y2x_numbering(bins_y{kk}(ll));
        
        voxel2panel_y(unk_id,1)=ind_panel; % front panel
        
        ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
        
        voxel2panel_y(unk_id,2)=ind_panel; % back panel
    end
    
end

clear map_from_y2x_numbering bins_y
num_panel_y=max(max(voxel2panel_y));

% Numbering panels along z direction
ind_panel=0; % counter for all panels along z
voxel2panel_z=zeros(num_nonair_cube,2); % 1 bottom, 2 up
for kk=1:length(bins_z)
    
    for ll=1:length(bins_z{kk})
        
        if (ll == 1) % beginning of connected elements, put a panel
            ind_panel=ind_panel+1;
        end
        
        unk_id=map_from_z2x_numbering(bins_z{kk}(ll));
        
        voxel2panel_z(unk_id,1)=ind_panel; % front panel
        
        ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
        
        voxel2panel_z(unk_id,2)=ind_panel; % back panel
    end
    
end
clear map_from_z2x_numbering bins_z

num_panel_z=max(max(voxel2panel_z));

voxel2panel_y = num_panel_x + voxel2panel_y;

voxel2panel_z = (num_panel_x+num_panel_y) + voxel2panel_z;

all_panels_ids=[voxel2panel_x; voxel2panel_y; voxel2panel_z];

clear voxel2panel_x voxel2panel_y voxel2panel_z

if(fl_profile == 1); disp(['Time for indexing panels ::: ', num2str(toc)]); end;

% The following is for visualization!!!
% retrieving the centers of non-air voxels
tic
dum=1;
xy_curr_ids_locs=zeros(3*num_nonair_cube,3);
for mm=1:N
    for ll=1:M
        for kk=1:L
            if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                    grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                xy_curr_ids_locs(dum,1:3)=coor_tmp;
                xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                xy_curr_ids_locs(2*num_nonair_cube+dum,1:3)=coor_tmp;
                dum=dum+1;
            end
        end
    end
end
if(fl_profile == 1); disp(['Time for retrieving centers of non-air voxels ::: ', num2str(toc)]);end;


% finding locations of panels enclosing x, y, and z currents
tic
all_panel_locs=zeros(3*num_nonair_cube,6);
for kk=1:num_nonair_cube
    %Jx currents
    all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
    all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
    %Jy currents
    ind=num_nonair_cube+kk;
    all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
    all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
    %Jz currents
    ind2=2*num_nonair_cube+kk;
    all_panel_locs(ind2,1:3)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)-dx*0.5];
    all_panel_locs(ind2,4:6)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)+dx*0.5];
end
clear xy_curr_ids_locs


% figure;
% plot3(all_panel_locs(:,1),all_panel_locs(:,2),all_panel_locs(:,3),'bo');
% hold on
% plot3(all_panel_locs(:,1),all_panel_locs(:,2),all_panel_locs(:,3),'r+');

if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end

% 8) Extract unique panel locations (boundary and inner)- and only boundary 
% panel locations and connectivity between them

% unique_panel_locs stores the locations and orientations of unique panels and to which
% each panel belongs to.
% (:,1:3) the locations
% (:,4) the orientations: 1-> x-aligned, 2->y-aligned, 3-> z-aligned; 
% (:,5) global id of panel 
% (:,6) conductor number


num_unique_panel=max(max(all_panels_ids));
unique_panel_locs=zeros(num_unique_panel,6);
for kk=1:size(all_panels_ids,1)
    
    cond_id=idxS_cond_dum(kk);
    
    unique_panel_locs(all_panels_ids(kk,1),1:3)=all_panel_locs(kk,1:3);
    
    if (all_panels_ids(kk,1) <= num_panel_x)
        unique_panel_locs(all_panels_ids(kk,1),4)=1; % x-aligned
    elseif (all_panels_ids(kk,1) <= (num_panel_x+num_panel_y))
        unique_panel_locs(all_panels_ids(kk,1),4)=2; % y-aligned
    else
        unique_panel_locs(all_panels_ids(kk,1),4)=3; % z-aligned
    end
    
    unique_panel_locs(all_panels_ids(kk,1),6)=cond_id;
    
    unique_panel_locs(all_panels_ids(kk,2),1:3)=all_panel_locs(kk,4:6);
    
    if (all_panels_ids(kk,2) <= num_panel_x)
        unique_panel_locs(all_panels_ids(kk,2),4)=1; % x-aligned
    elseif (all_panels_ids(kk,2) <= (num_panel_x+num_panel_y))
        unique_panel_locs(all_panels_ids(kk,2),4)=2; % y-aligned
    else
        unique_panel_locs(all_panels_ids(kk,2),4)=3; % z-aligned
    end
    
    unique_panel_locs(all_panels_ids(kk,2),6)=cond_id;

end
%clear all_panel_locs

% find out whether shared panel or boundary panel
pnl_shared_or_not=zeros(num_unique_panel,1);
for kk=1:size(all_panels_ids,1)
    pnl_shared_or_not(all_panels_ids(kk,1))=pnl_shared_or_not(all_panels_ids(kk,1))+1;
    pnl_shared_or_not(all_panels_ids(kk,2))=pnl_shared_or_not(all_panels_ids(kk,2))+1;
end

%pnl_shared_id=find(pnl_shared_or_not==2);
pnl_bndry_id=find(pnl_shared_or_not==1);
unique_panel_bndry_locs=zeros(length(pnl_bndry_id),6);

for kk=1:length(pnl_bndry_id)
   unique_panel_bndry_locs(kk,1:6)=unique_panel_locs(pnl_bndry_id(kk),1:6);
end

% 9) Finding the ids of panels in global grid 

bbox_origin=[min(x_bnd) min(y_bnd) min(z_bnd)];
bbox_elem(1)=x_bnd(2)-x_bnd(1);
bbox_elem(2)=y_bnd(2)-y_bnd(1);
bbox_elem(3)=z_bnd(2)-z_bnd(1);
for kk=1:3
    nelem_xyz(kk)=floor(bbox_elem(kk)/dx);
end

% generate ids of x_directed panels
tensor_ijk_xdir_pnl=zeros(nelem_xyz(1)+1,nelem_xyz(2),nelem_xyz(3));
dum_cnt=1;
for mm=1:nelem_xyz(3)% z-variation
    for ll=1:nelem_xyz(2) % y-variation
        for kk=1:nelem_xyz(1)+1 % x-variation
            tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end

% generate ids of  y_directed panels
tensor_ijk_ydir_pnl=zeros(nelem_xyz(1),nelem_xyz(2)+1,nelem_xyz(3));
for kk=1:nelem_xyz(1) % x-variation
    for mm=1:nelem_xyz(3) % z-variation
        for ll=1:nelem_xyz(2)+1 % y-variation
            tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end

% generate ids of z_directed panels
tensor_ijk_zdir_pnl=zeros(nelem_xyz(1),nelem_xyz(2),nelem_xyz(3)+1);
for ll=1:nelem_xyz(2) % y-variation
    for kk=1:nelem_xyz(1) % x-variation
        for mm=1:nelem_xyz(3)+1 % z-variation
            tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end

% assign these ids to the unique_panels_loc
org_x_panel=[bbox_origin(1) bbox_origin(2)+dx/2 bbox_origin(3)+dx/2];
org_y_panel=[bbox_origin(1)+dx/2 bbox_origin(2) bbox_origin(3)+dx/2];
org_z_panel=[bbox_origin(1)+dx/2 bbox_origin(2)+dx/2 bbox_origin(3)];

% for unique panels
for kk=1:size(unique_panel_locs,1)
   if (unique_panel_locs(kk,4) == 1) % x-directed panel
       tmp_ind=round((unique_panel_locs(kk,1:3)-org_x_panel)/dx)+1;
       unique_panel_locs(kk,5)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   elseif (unique_panel_locs(kk,4) == 2) % y-directed panel
       tmp_ind=round((unique_panel_locs(kk,1:3)-org_y_panel)/dx)+1;
       unique_panel_locs(kk,5)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   elseif (unique_panel_locs(kk,4) == 3) % z-directed panel
       tmp_ind=round((unique_panel_locs(kk,1:3)-org_z_panel)/dx)+1;
       unique_panel_locs(kk,5)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   end
end

% for unique boundary panels
for kk=1:size(unique_panel_bndry_locs,1)
   if (unique_panel_bndry_locs(kk,4) == 1) % x-directed panel
       tmp_ind=round((unique_panel_bndry_locs(kk,1:3)-org_x_panel)/dx)+1;
       unique_panel_bndry_locs(kk,5)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   elseif (unique_panel_bndry_locs(kk,4) == 2) % y-directed panel
       tmp_ind=round((unique_panel_bndry_locs(kk,1:3)-org_y_panel)/dx)+1;
       unique_panel_bndry_locs(kk,5)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   elseif (unique_panel_bndry_locs(kk,4) == 3) % z-directed panel
       tmp_ind=round((unique_panel_bndry_locs(kk,1:3)-org_z_panel)/dx)+1;
       unique_panel_bndry_locs(kk,5)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
   end
end


% 10) Visualization: check whether the data generated correctly

% 10a) visualize graphs along x, y, or z directions
fl_vis_graphs=0;
if(fl_vis_graphs == 1)
    figure; plot(G_x);title('Graph for panels along x direction')
    figure; plot(G_y);title('Graph for panels along y direction')
    figure; plot(G_z);title('Graph for panels along z direction')
end

% 10b) visualize the numbers of boxes for Jx, Jy, and Jz currents
fl_vis_numbering=0;
if (fl_vis_numbering == 1)
    figure;
    set(gca,'FontSize',24);
    xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');set(gca,'FontSize',24);
    hold on
    for mm=1:N
        for ll=1:M
            for kk=1:L
                coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,1))); % x numbering
                %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,2))); % y numbering
                %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,3))); % z numbering
                set(h,'FontSize',24)
            end
        end
    end
    xlabel('x');ylabel('y'); set(gca,'FontSize',24);
    grid on
    axis tight
    grid on
    view(2)
end

% 10c) check whether the panel locations are correct
if (fl_check_pnls > 0)
    
    % select what you'll plot - all panels or only boundary panels - one of
    % the following:
    if (fl_check_pnls  == 2)
        tmp_mat=unique_panel_bndry_locs;
    elseif (fl_check_pnls  == 1)
        tmp_mat=unique_panel_locs;        
    end
    
    figure
    for kk=1:size(tmp_mat,1)
        verts=[];
        if (tmp_mat(kk,4) == 1) % x-directed
            verts=[tmp_mat(kk,1) tmp_mat(kk,2)-dx/2 tmp_mat(kk,3)-dx/2; ...
                tmp_mat(kk,1) tmp_mat(kk,2)+dx/2 tmp_mat(kk,3)-dx/2; ...
                tmp_mat(kk,1) tmp_mat(kk,2)+dx/2 tmp_mat(kk,3)+dx/2; ...
                tmp_mat(kk,1) tmp_mat(kk,2)-dx/2 tmp_mat(kk,3)+dx/2;];
        elseif (tmp_mat(kk,4) == 2) % y-directed
            verts=[tmp_mat(kk,1)-dx/2 tmp_mat(kk,2) tmp_mat(kk,3)-dx/2; ...
                tmp_mat(kk,1)+dx/2 tmp_mat(kk,2) tmp_mat(kk,3)-dx/2; ...
                tmp_mat(kk,1)+dx/2 tmp_mat(kk,2) tmp_mat(kk,3)+dx/2; ...
                tmp_mat(kk,1)-dx/2 tmp_mat(kk,2) tmp_mat(kk,3)+dx/2;];
        elseif (tmp_mat(kk,4) == 3) % z-directed
            verts=[tmp_mat(kk,1)-dx/2 tmp_mat(kk,2)-dx/2 tmp_mat(kk,3); ...
                tmp_mat(kk,1)+dx/2 tmp_mat(kk,2)-dx/2 tmp_mat(kk,3); ...
                tmp_mat(kk,1)+dx/2 tmp_mat(kk,2)+dx/2 tmp_mat(kk,3); ...
                tmp_mat(kk,1)-dx/2 tmp_mat(kk,2)+dx/2 tmp_mat(kk,3);];
        end
        patch(verts(:,1),verts(:,2),verts(:,3),[1 2 4 3],'EdgeColor','blue','FaceColor','none');
        hold on
        %plot3(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),'r+')
        %hold on
        %text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(kk,5))
        text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(tmp_mat(kk,5)))
        %text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(tmp_mat(kk,6)))
        hold on
    end
    xlabel('x');ylabel('y');zlabel('z');view(-30,30);set(gca,'FontSize',16); axis tight;
end


disp('Done... Extracting unique panels of a grid')
disp('-----------------------------------------------------')




