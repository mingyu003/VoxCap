clc; close all; clear all; format long e;
% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------

pre_define_the_path_for_folders

% -------------------------------------------------------------------------
%                  Inputs for Simulation
% -------------------------------------------------------------------------
er = 0;  % epsilon_r of conductors
epsa=[5 2.6 5 2.6 5 2.6 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7 3.7];% dielectric permittivity, for conductor case, set it to null []
epsb=1; % background permittivity
eps0=8.854187817e-12;
se=5.8e7; % conductivity of conductors
inner_it = 50; outer_it = 20; tol=1e-4; % iterative solver inputs
Res = 10 ; % voxel size (deltax)
tole=1e-8;%tolerance for Tucker decompression
fl_check_domain=0; % set to 1 for only plotting the structure
fl_check_geo=0; % set to 1 for only plotting the domain
fl_check_domain_pnls=0; % set to 1 for  boundary panels 
fl_save_data = 1; % set 1 for saving data for post-processing later on
fl_write_fastcap_file = 0; % set to 1 for outputing the panel info for fastcap sim
fl_Tucker_decomp=0; % set 1 to compress tensor with Tucker
fl_precond=2; %set 1 to implement preconditioner, 2 to preconditioner_ver2_n
num_vox_in_blk=10;% # of voxels in each preconditioner block. 
fl_filling_or_retrieval=0; % set to 1 for filling the circulant tensor and 0 for retrieval from prestored Toeplitz tensor
% -------------------------------------------------------------------------
%                  Inputs for the Structure
% -------------------------------------------------------------------------
% We only need centers (Cnt), dimensions (Dims), orientations (Orients),
% inner/outer permittivities (Eps_inout) of conductors/dielectrics
% at the end of this part. Assumption: inner permittivity of conductor is
% zero & reference points to find the directions of panel normals are
% center points of dielectrics.

% predefine the #s of conductors and dielectrics
num_cond_x=15;
num_cond_z=3;
num_cond=num_cond_x*num_cond_z;

num_diel=length(epsa);
% cross-sections and lengths of buses
width_bus=70;
height_bus=140;
len_bus=2030;

dist=2*width_bus; % between centers of 2 buses
dist_z=height_bus+80+40;

if num_diel>0
    buff_btw_cond_diel=width_bus; % buffer to put dielectric around the conductors
    shft_vect=[buff_btw_cond_diel buff_btw_cond_diel buff_btw_cond_diel];
else
    shft_vect=[0 0 0];
end

%%%%%%%%%% cross buses section
% inputs for conductor
cen_cond=zeros(num_cond,3);
counter=1;
for kk=1:num_cond_z
    for ll=1:num_cond_x
        if kk==1
            cen_cond(counter,1:3)=[len_bus/2+shft_vect(1) width_bus/2+(dist*(ll-1))+shft_vect(2) 40+height_bus/2];
        elseif kk==2
            cen_cond(counter,1:3)=[(width_bus+0)/2+((dist-0)*(ll-1))+shft_vect(1) len_bus/2+shft_vect(2) 300+height_bus/2];
        elseif kk==3
            cen_cond(counter,1:3)=[len_bus/2+shft_vect(1) width_bus/2+(dist*(ll-1))+shft_vect(2) 560+height_bus/2];
        end
%         cen_cond(2*counter-1,1:3)=[width_bus/2+(dist*(ll-1)) len_bus/2 height_bus/2+(dist*(kk-1))] + num_diel*shft_vect;
%         cen_cond(2*counter,1:3)=[width_bus/2+(dist*(ll-1)) len_bus/2 height_bus/2+3+(dist*(kk-1))] + num_diel*shft_vect;
        counter=counter+1;
    end
end


counter=1;
Cnt_cond=[cen_cond];Ref_pnts_cond=[];Eps_inout_cond=[];Dims_cond=[];Orients_cond=[];
for kk=1:num_cond_z
    for ll=1:num_cond_x
        %         Cnt_cond = [Cnt_cond; cen_cond(counter,1:3);cen_cond(counter+1,1:3)];
%         Ref_pnts_cond = [Ref_pnts_cond; cen_cond(counter,1:3);cen_cond(counter+1,1:3)];
        if (mod(kk,2) == 1)
            tmp = [len_bus width_bus   height_bus];
        else
            tmp = [len_bus (width_bus+0)   height_bus];
        end
        Dims_cond = [Dims_cond; tmp];
        %         Orients_cond=[Orients_cond;'x';'x'];
        if (mod(kk,2) == 1)
            Orients_cond=[Orients_cond;'x'];
        else
            Orients_cond=[Orients_cond;'y'];
        end
        counter=counter+1;
    end
end

%%%%%%%%%%%%% dielectric section
% inputs for dielectric
if num_diel>0
    Dims_diel=zeros(num_diel,3);Cnt_diel=zeros(num_diel,3);Eps_inout_diel=zeros(num_diel,2);Orients_diel=[];
    % diel 1: layer 2
    Dims_diel(1,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(1,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(1,3)=40; % height of dielectric-ii
    Cnt_diel(1,:)=[Dims_diel(1,1)/2  Dims_diel(1,2)/2  Dims_diel(1,3)/2];
    Orients_diel=[Orients_diel;'y';];
    % diel 3: layer 4
    Dims_diel(3,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(3,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(3,3)=40; % height of dielectric-ii
    Cnt_diel(3,:)=[Dims_diel(3,1)/2  Dims_diel(3,2)/2  Dims_diel(3,3)/2+260];
    Orients_diel=[Orients_diel;'y';];
    % diel 5: layer 6
    Dims_diel(5,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(5,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(5,3)=40; % height of dielectric-ii
    Cnt_diel(5,:)=[Dims_diel(5,1)/2  Dims_diel(5,2)/2  Dims_diel(5,3)/2+520];
    Orients_diel=[Orients_diel;'y';];
    
    
    % diel 2: layer 3
    Dims_diel(2,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(2,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(2,3)=220; % height of dielectric-ii
    Cnt_diel(2,:)=[Dims_diel(2,1)/2  Dims_diel(2,2)/2  Dims_diel(2,3)/2+40];
    Orients_diel=[Orients_diel;'y';];
    % diel 4: layer 5
    Dims_diel(4,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(4,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(4,3)=220; % height of dielectric-ii
    Cnt_diel(4,:)=[Dims_diel(4,1)/2  Dims_diel(4,2)/2  Dims_diel(4,3)/2+300];
    Orients_diel=[Orients_diel;'y';];
    % diel 6: layer 7
    Dims_diel(6,1)=len_bus+2*shft_vect(1); % length of dielectric_ii
    Dims_diel(6,2)=len_bus+2*shft_vect(2); % width of dielectric-ii
    Dims_diel(6,3)=220; % height of dielectric-ii
    Cnt_diel(6,:)=[Dims_diel(6,1)/2  Dims_diel(6,2)/2  Dims_diel(6,3)/2+560];
    Orients_diel=[Orients_diel;'y';];
    
    for ii=7:21
        Dims_diel(ii,1)=len_bus+20; % length of dielectric_ii
        Dims_diel(ii,2)=width_bus+20+0; % width of dielectric-ii
        Dims_diel(ii,3)=height_bus+20; % height of dielectric-ii
        Cnt_diel(ii,:)=cen_cond(ii+9,:);
        Orients_diel=[Orients_diel;'y';];
    end
    
end

if num_diel>0
    Ref_pnts_diel=Cnt_diel;
    Cnt=[Cnt_cond;Cnt_diel];
    Ref_pnts=[Ref_pnts_cond;Ref_pnts_diel];
    Eps_inout=[2.6 5];
else
    Cnt=Cnt_cond;
    Ref_pnts=Cnt_cond;
    Eps_inout=[2.6 5];
end
if num_diel>0
    Dims=[Dims_cond;Dims_diel];
    Orients=[Orients_cond;Orients_diel];
else
    Dims=Dims_cond;
    Orients=Orients_cond;
end
% -------------------------------------------------------------------------
%                  Input for Computational Domain
% -------------------------------------------------------------------------
% At the end of this part, we only need bbox_min(3) and bbox_max(3) vectors
% define computational domain or bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len_bus+dist len_bus+dist 780]+1e-8 ; % max coordinates of bbox

% bbox_max=[width_bus width_bus height_bus] ; % for conductor case
% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------
pre_print_out_inputs_generate_consts
tini_pre = tic;
% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);
% assign constitutive parameters
[idx,grid_intcon] = mingyu_distinguishvoxelsofstructure(Res,r,Cnt,Dims,Orients);
% [idx,grid_intcon] = sphere_distinguishvoxelsofstructure(r,Res,Cnt,Dims,fl_check_geo);
% ------------------------------------------------------------------------
%                  Obtain Panel Coordinates and IDs
% -------------------------------------------------------------------------
[L,M,N,~] = size(r); % domain size
dx = Res; % voxel size

disp('Generating panels of voxelized computational domain ...')
grid_tmp = ones(L,M,N); idx_dum = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp;
idx_dum_cond_no=ones(length(idx_dum),1);

[comp_dom_panels]=new_computational_panels(dx,r);
clear r

disp('Generating panels of voxelized individual bodies ...')

geom_bndry_panels=cell(size(grid_intcon,1),1);
% geom_bndry_panels=cell(size(grid_intcon,1),1);
for kk=1:size(grid_intcon,1)
    idx_cd_no(1:length(idx{kk}))=kk;
    [geom_bndry_panels{kk}]=new_geom_panels(dx,grid_intcon{kk},idx{kk},idx_cd_no,num_vox_in_blk);
    clear idx_cd_no
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear grid_intcon idx;
clear grid_intcon

fl_check_geo_pnls=0;
if fl_check_geo_pnls>0
    if fl_check_geo_pnls==1
    plot_panel=cell2mat(geom_bndry_panels); 
%     else
%     plot_panel=cell2mat(geom_panels);
    end
    plot_panels(dx,plot_panel);
    clear plot_panel;
    return
end
if fl_check_domain_pnls==1
    plot_panels(dx,comp_dom_panels);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                Clear overlapped surfaces                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:15
    geom_bndry_panels{num_cond+6+ii}(find(geom_bndry_panels{num_cond+6+ii}(:,3)<300),:)=[];
end

% z_bottom_diel_2=find(geom_bndry_panels{32}(:,3)==0.04);
geom_bndry_panels{num_cond+2}(find(geom_bndry_panels{num_cond+2}(:,3)==40),:)=[];
% z_bottom_diel_3=find(geom_bndry_panels{33}(:,3)==0.26);
geom_bndry_panels{num_cond+3}(find(geom_bndry_panels{num_cond+3}(:,3)==260),:)=[];
% z_bottom_diel_4=find(geom_bndry_panels{34}(:,3)==0.3);
geom_bndry_panels{num_cond+4}(find(geom_bndry_panels{num_cond+4}(:,3)==300),:)=[];
geom_bndry_panels{num_cond+5}(find(geom_bndry_panels{num_cond+5}(:,3)==520),:)=[];
geom_bndry_panels{num_cond+6}(find(geom_bndry_panels{num_cond+6}(:,3)==560),:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        Define the 7 8 and 9th columns in geom_bndry_panels   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:num_cond_x
    geom_bndry_panels{ii}(:,7)=0;geom_bndry_panels{ii}(:,8)=0;geom_bndry_panels{ii}(:,9)=2.6;
end
for ii=num_cond_x+1:2*num_cond_x
    geom_bndry_panels{ii}(:,7)=0;geom_bndry_panels{ii}(:,8)=0;geom_bndry_panels{ii}(:,9)=3.7;
end
for ii=2*num_cond_x+1:3*num_cond_x
    geom_bndry_panels{ii}(:,7)=0;geom_bndry_panels{ii}(:,8)=0;geom_bndry_panels{ii}(:,9)=2.6;
end

geom_bndry_panels{num_cond+1}(:,7)=1;geom_bndry_panels{num_cond+1}(:,8)=5.0;geom_bndry_panels{num_cond+1}(:,9)=1;
geom_bndry_panels{num_cond+2}(:,7)=2;geom_bndry_panels{num_cond+2}(:,8)=2.6;geom_bndry_panels{num_cond+2}(:,9)=1;
geom_bndry_panels{num_cond+3}(:,7)=3;geom_bndry_panels{num_cond+3}(:,8)=5.0;geom_bndry_panels{num_cond+3}(:,9)=1;
geom_bndry_panels{num_cond+4}(:,7)=4;geom_bndry_panels{num_cond+4}(:,8)=2.6;geom_bndry_panels{num_cond+4}(:,9)=1;
geom_bndry_panels{num_cond+5}(:,7)=5;geom_bndry_panels{num_cond+5}(:,8)=5.0;geom_bndry_panels{num_cond+5}(:,9)=1;
geom_bndry_panels{num_cond+6}(:,7)=6;geom_bndry_panels{num_cond+6}(:,8)=2.6;geom_bndry_panels{num_cond+6}(:,9)=1;

for ii=num_cond+6+1:num_cond+6+15 %%add for merge M!M2 && coat M2
    geom_bndry_panels{ii}(:,7)=ii-45;geom_bndry_panels{ii}(:,8)=3.7;geom_bndry_panels{ii}(:,9)=2.6;%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%% interface diel1 && diel2
z_up_diel1=find(geom_bndry_panels{num_cond+1}(:,3)==40);
z_up_diel1_cond_interface=cell(num_cond_x,1);
for ii=1:num_cond_x
    z_up_diel1_cond_interface{ii}=find(geom_bndry_panels{num_cond+1}(:,3)==40 & geom_bndry_panels{num_cond+1}(:,1)>=70 & geom_bndry_panels{num_cond+1}(:,1)<=2100 & ...
       geom_bndry_panels{num_cond+1}(:,2)>=70+(ii-1)*140 & geom_bndry_panels{num_cond+1}(:,2)<=140+(ii-1)*140);
end
z_up_diel1_diel2_interface=setdiff(z_up_diel1,cell2mat(z_up_diel1_cond_interface));
geom_bndry_panels{num_cond+1}(z_up_diel1_diel2_interface,7)=num_diel+1;%%%%%%%%%%%%%%%%
geom_bndry_panels{num_cond+1}(z_up_diel1_diel2_interface,8)=5;
geom_bndry_panels{num_cond+1}(z_up_diel1_diel2_interface,9)=2.6;
geom_bndry_panels{num_cond+1}(cell2mat(z_up_diel1_cond_interface),:)=[];

%%%%%% interface diel2 && diel3
z_up_diel2=find(geom_bndry_panels{num_cond+2}(:,3)==260);
geom_bndry_panels{num_cond+2}(z_up_diel2,7)=num_diel+2; %%%%%%%%%%%%%%
geom_bndry_panels{num_cond+2}(z_up_diel2,8)=2.6;
% geom_bndry_panels{20}(z_up_diel2,9)=2.6; %add for merge M1M2
geom_bndry_panels{num_cond+2}(z_up_diel2,9)=5;

%%%%%% interface diel3 && cond 7-12 , diel3 && diel4
z_up_diel3=find(geom_bndry_panels{num_cond+3}(:,3)==300);
z_up_diel3_cond_interface=cell(num_cond_x,1);
for ii=1:num_cond_x
        z_up_diel3_cond_interface{ii}=find(geom_bndry_panels{num_cond+3}(:,3)==300 & geom_bndry_panels{num_cond+3}(:,2)>=70 & geom_bndry_panels{num_cond+3}(:,2)<=2100 & ...
           geom_bndry_panels{num_cond+3}(:,1)>=70+(ii-1)*140 & geom_bndry_panels{num_cond+3}(:,1)<=140+(ii-1)*140);
%     z_up_diel3_cond_interface{ii}=find(geom_bndry_panels{num_cond+3}(:,3)==300 & geom_bndry_panels{num_cond+3}(:,2)>=70 & geom_bndry_panels{num_cond+3}(:,2)<=1400 & ...
%         geom_bndry_panels{num_cond+3}(:,1)>=70+(ii-1)*130 & geom_bndry_panels{num_cond+3}(:,1)<=170+(ii-1)*130);
end
z_up_diel3_diel4_interface=setdiff(z_up_diel3,cell2mat(z_up_diel3_cond_interface));
geom_bndry_panels{num_cond+3}(z_up_diel3_diel4_interface,7)=num_diel+3;
% geom_bndry_panels{21}(z_up_diel3_diel4_interface,8)=2.6; %add for merge M1M2
geom_bndry_panels{num_cond+3}(z_up_diel3_diel4_interface,8)=5;
geom_bndry_panels{num_cond+3}(z_up_diel3_diel4_interface,9)=2.6; 
geom_bndry_panels{num_cond+3}(cell2mat(z_up_diel3_cond_interface),:)=[];

%%%%%% interface diel4 && diel5
z_up_diel4=find(geom_bndry_panels{num_cond+4}(:,3)==520);
geom_bndry_panels{num_cond+4}(z_up_diel4,7)=num_diel+4;%%%%%%%%%%%%%%%%%%%%
geom_bndry_panels{num_cond+4}(z_up_diel4,8)=2.6;
geom_bndry_panels{num_cond+4}(z_up_diel4,9)=5;

%%%%%% interface diel5 && diel6
z_up_diel5=find(geom_bndry_panels{num_cond+5}(:,3)==560);
z_up_diel5_cond_interface=cell(num_cond_x,1);
for ii=1:num_cond_x
    z_up_diel5_cond_interface{ii}=find(geom_bndry_panels{num_cond+5}(:,3)==560 & geom_bndry_panels{num_cond+5}(:,1)>=70 & geom_bndry_panels{num_cond+5}(:,1)<=2100 & ...
       geom_bndry_panels{num_cond+5}(:,2)>=70+(ii-1)*140 & geom_bndry_panels{num_cond+5}(:,2)<=140+(ii-1)*140);
end
z_up_diel5_diel6_interface=setdiff(z_up_diel5,cell2mat(z_up_diel5_cond_interface));
geom_bndry_panels{num_cond+5}(z_up_diel5_diel6_interface,7)=num_diel+5;
geom_bndry_panels{num_cond+5}(z_up_diel5_diel6_interface,8)=5;
geom_bndry_panels{num_cond+5}(z_up_diel5_diel6_interface,9)=2.6;
geom_bndry_panels{num_cond+5}(cell2mat(z_up_diel5_cond_interface),:)=[];

%-----------------------------------------------------------
fl_plot_pnl_orients=0;
geom_bndry_panels=mingyu_pre_correct_pnl_orients(num_cond,num_diel,geom_bndry_panels,Eps_inout,Cnt,dx,fl_plot_pnl_orients);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% interface information  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_up_diel1_diel2_interface=find(geom_bndry_panels{num_cond+1}(:,3)==40);
% geom_bndry_panels{31}=geom_bndry_panels{19}(z_up_diel1_diel2_interface,:); %interface  %add for merge M1M2
% geom_bndry_panels{31}(:,6)=num_cond+num_diel+1;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+num_diel+1}=geom_bndry_panels{num_cond+1}(z_up_diel1_diel2_interface,:); %interface 
geom_bndry_panels{num_cond+num_diel+1}(:,6)=num_cond+num_diel+1;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+1}(z_up_diel1_diel2_interface,:)=[];%delete overlapped surface

geom_bndry_panels{num_cond+num_diel+2}=geom_bndry_panels{num_cond+2}(z_up_diel2,:); %interface
geom_bndry_panels{num_cond+num_diel+2}(:,6)=num_cond+num_diel+2;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+2}(z_up_diel2,:)=[];%delete overlapped surface

z_up_diel3_diel4_interface=find(geom_bndry_panels{num_cond+3}(:,3)==300);
geom_bndry_panels{num_cond+num_diel+3}=geom_bndry_panels{num_cond+3}(z_up_diel3_diel4_interface,:); %interface
geom_bndry_panels{num_cond+num_diel+3}(:,6)=num_cond+num_diel+3;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+3}(z_up_diel3_diel4_interface,:)=[];%delete overlapped surface

geom_bndry_panels{num_cond+num_diel+4}=geom_bndry_panels{num_cond+4}(z_up_diel4,:); %interface
geom_bndry_panels{num_cond+num_diel+4}(:,6)=num_cond+num_diel+4;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+4}(z_up_diel4,:)=[];%delete overlapped surface
% 
z_up_diel5_diel6_interface=find(geom_bndry_panels{num_cond+5}(:,3)==560);
geom_bndry_panels{num_cond+num_diel+5}=geom_bndry_panels{num_cond+5}(z_up_diel5_diel6_interface,:); %interface
geom_bndry_panels{num_cond+num_diel+5}(:,6)=num_cond+num_diel+5;% multi-interfaces, change 1 to ii
geom_bndry_panels{num_cond+5}(z_up_diel5_diel6_interface,:)=[];%delete overlapped surface

% geom_bndry_panels{32}=geom_bndry_panels{22}(z_up_diel4,:); %interface
% geom_bndry_panels{32}(:,6)=num_cond+num_diel+2;% multi-interfaces, change 1 to ii
% geom_bndry_panels{32}(:,7)=14;% multi-interfaces, change 1 to ii
% geom_bndry_panels{22}(z_up_diel4,:)=[];%delete overlapped surface
% 
% z_up_diel5_diel6_interface=find(geom_bndry_panels{23}(:,3)==560);
% geom_bndry_panels{33}=geom_bndry_panels{23}(z_up_diel5_diel6_interface,:); %interface
% geom_bndry_panels{33}(:,6)=num_cond+num_diel+3;% multi-interfaces, change 1 to ii
% geom_bndry_panels{33}(:,7)=15;% multi-interfaces, change 1 to ii
% geom_bndry_panels{23}(z_up_diel5_diel6_interface,:)=[];%delete overlapped surface
geom_bndry_panels3=geom_bndry_panels;
%  return


num_interface=5;
% num_interface=3;

ss=zeros(num_diel+num_interface,1);
for ii=1:num_diel
    epsa=geom_bndry_panels{num_cond+ii}(1,8);
    epsb=geom_bndry_panels{num_cond+ii}(1,9);
    ss(ii)=(dx^2)*((epsa+epsb)/((epsa-epsb)*2*eps0));
end
for ii=1:num_interface
    epsa=geom_bndry_panels{num_cond+num_diel+ii}(1,8);
    epsb=geom_bndry_panels{num_cond+num_diel+ii}(1,9);
    ss(num_diel+ii)=(dx^2)*((epsa+epsb)/((epsa-epsb)*2*eps0));
end


if (num_diel == 0)
    fl_st_diel=1;
else
    fl_st_diel=num_diel+num_cond+num_interface;
end

num_panels=cell(fl_st_diel,1);
ids_panels=cell(fl_st_diel,1);
num_tot_panels=0;
geom_bndry_panels_cond=cell(num_cond,1);
for ii=1:num_cond
    geom_bndry_panels_cond{ii}=cell2mat(geom_bndry_panels(ii));
    geom_bndry_panels_cond{ii}=sortrows(geom_bndry_panels_cond{ii},5);
end

for kk=1:num_cond
    [num_panels{kk},ids_panels{kk}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels_cond{kk});
    num_tot_panels=num_tot_panels+num_panels{kk};
end

if num_diel>0
    for kk=1:num_diel
        [num_panels{kk+num_cond},ids_panels{kk+num_cond}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels{kk+num_cond});
        num_tot_panels=num_tot_panels+num_panels{kk+1};
    end
end

for kk=1:num_interface
    [num_panels{kk+num_cond+num_diel},ids_panels{kk+num_cond+num_diel}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels{kk+num_cond+num_diel});
    num_tot_panels=num_tot_panels+num_panels{kk+1};
end

clear comp_dom_panels 

if (num_diel == 0)
    ids_panels{2}=cell(6,1);
end

geom_bndry_panels2=cell(fl_st_diel,1);
for ii=1:num_cond
    geom_bndry_panels2{ii}=geom_bndry_panels_cond{ii};
end
clear geom_bndry_panels_cond

if(num_diel > 0)% just add at the end
    for ii=1:num_diel
        geom_bndry_panels2{ii+num_cond} = geom_bndry_panels{num_cond+ii};
    end
end
for ii=1:num_interface
    geom_bndry_panels2{ii+num_cond+num_diel} = geom_bndry_panels{num_cond+num_diel+ii};
end
clear  geom_bndry_panels

geom_bndry_panels=zeros(num_tot_panels,10);
st_ind=1; end_ind=0;
for kk=1:fl_st_diel
    end_ind=end_ind+num_panels{kk};
    geom_bndry_panels(st_ind:end_ind,1:10)=geom_bndry_panels2{kk}(:,1:10);
    st_ind=end_ind+1;
end
clear  geom_bndry_panels2

diel_pnl_id_locs=cell(num_diel+num_interface,1);diel_pnl_ids=cell(num_diel+num_interface,1);
for ii=1:num_diel
    diel_pnl_id_locs{ii}=cell(6,1); diel_pnl_ids{ii}=cell(6,1);
    if (num_diel > 0) % currently num_diel could be max 1
        % information for +/- direction directed panels from combined ids_panels
        diel_pnl_id_locs{ii}{1}=find(ids_panels{num_cond+ii}{4}(:)>0);
        diel_pnl_ids{ii}{1}=ids_panels{num_cond+ii}{1}(diel_pnl_id_locs{ii}{1});
        diel_pnl_id_locs{ii}{4}=find(ids_panels{num_cond+ii}{4}(:)<0);
        diel_pnl_ids{ii}{4}=ids_panels{num_cond+ii}{1}(diel_pnl_id_locs{ii}{4});
        
        diel_pnl_id_locs{ii}{2}=find(ids_panels{num_cond+ii}{5}(:)>0);
        diel_pnl_ids{ii}{2}=ids_panels{num_cond+ii}{2}(diel_pnl_id_locs{ii}{2});
        diel_pnl_id_locs{ii}{5}=find(ids_panels{num_cond+ii}{5}(:)<0);
        diel_pnl_ids{ii}{5}=ids_panels{num_cond+ii}{2}(diel_pnl_id_locs{ii}{5});
        
        diel_pnl_id_locs{ii}{3}=find(ids_panels{num_cond+ii}{6}(:)>0);
        diel_pnl_ids{ii}{3}=ids_panels{num_cond+ii}{3}(diel_pnl_id_locs{ii}{3});
        diel_pnl_id_locs{ii}{6}=find(ids_panels{num_cond+ii}{6}(:)<0);
        diel_pnl_ids{ii}{6}=ids_panels{num_cond+ii}{3}(diel_pnl_id_locs{ii}{6});
        
        % diel_pnl_id_locs stores the locations of elements of diel_pnl_ids which
        % are positive/negative x, y, and z directed panels.
        % Its entries {1} +x-directed, {2} +y-directed, {3} +z-directed, {4}
        % -x-directed, {5} -y directed, and {6} -z-directed
    end
end
for ii=1:num_interface
%     ids_panels{num_diel+num_cond+ii}{6}(:)=-1;
    diel_pnl_id_locs{num_diel+ii}=cell(6,1); diel_pnl_ids{num_diel+ii}=cell(6,1);
%     if (num_diel > 0) % currently num_diel could be max 1
        % information for +/- direction directed panels from combined ids_panels
%         diel_pnl_id_locs{ii}{1}=find(ids_panels{num_cond+ii}{4}(:)>0);
%         diel_pnl_ids{ii}{1}=ids_panels{num_cond+ii}{1}(diel_pnl_id_locs{ii}{1});
%         diel_pnl_id_locs{ii}{4}=find(ids_panels{num_cond+ii}{4}(:)<0);
%         diel_pnl_ids{ii}{4}=ids_panels{num_cond+ii}{1}(diel_pnl_id_locs{ii}{4});
%         
%         diel_pnl_id_locs{ii}{2}=find(ids_panels{num_cond+ii}{5}(:)>0);
%         diel_pnl_ids{ii}{2}=ids_panels{num_cond+ii}{2}(diel_pnl_id_locs{ii}{2});
%         diel_pnl_id_locs{ii}{5}=find(ids_panels{num_cond+ii}{5}(:)<0);
%         diel_pnl_ids{ii}{5}=ids_panels{num_cond+ii}{2}(diel_pnl_id_locs{ii}{5});
        
        diel_pnl_id_locs{num_diel+ii}{3}=find(ids_panels{num_diel+num_cond+ii}{6}(:)>0);
        diel_pnl_ids{num_diel+ii}{3}=ids_panels{num_diel+num_cond+ii}{3}(diel_pnl_id_locs{num_diel+ii}{3});
        diel_pnl_id_locs{num_diel+ii}{6}=find(ids_panels{num_diel+num_cond+ii}{6}(:)<0);
        diel_pnl_ids{num_diel+ii}{6}=ids_panels{num_diel+num_cond+ii}{3}(diel_pnl_id_locs{num_diel+ii}{6});
        
        % diel_pnl_id_locs stores the locations of elements of diel_pnl_ids which
        % are positive/negative x, y, and z directed panels.
        % Its entries {1} +x-directed, {2} +y-directed, {3} +z-directed, {4}
        % -x-directed, {5} -y directed, and {6} -z-directed
%     end
end

% Attention: the following lines were added to distinguish x-, y-, z- aligned
% conductor and dielectric panels

inds_glob=zeros(3*(num_cond+num_diel+num_interface),2);
x_aligned_cond_inds=zeros(num_cond,1);
y_aligned_cond_inds=zeros(num_cond,1);
z_aligned_cond_inds=zeros(num_cond,1);
x_aligned_diel_inds=zeros(num_diel,1);
y_aligned_diel_inds=zeros(num_diel,1);
z_aligned_diel_inds=zeros(num_diel,1);
x_aligned_interface_inds=zeros(num_interface,1);
y_aligned_interface_inds=zeros(num_interface,1);
z_aligned_interface_inds=zeros(num_interface,1);

for ii=1:num_cond
    x_aligned_cond_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii));
    y_aligned_cond_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii));
    z_aligned_cond_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii));
end
for ii=1:num_diel
    x_aligned_diel_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == ii) );
    y_aligned_diel_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == ii) );
    z_aligned_diel_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == ii) );
end
for ii=1:num_interface
    x_aligned_interface_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == num_diel+ii) );
    y_aligned_interface_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == num_diel+ii) );
    z_aligned_interface_inds(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == num_diel+ii) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%   x-aligned  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inds_glob(1,1)=1;
inds_glob(1,2)=x_aligned_cond_inds(1);
for ii=2:num_cond
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+x_aligned_cond_inds(ii);
end
inds_glob(num_cond+1,1)=1+inds_glob(num_cond,2);
inds_glob(num_cond+1,2)=inds_glob(num_cond,2)+x_aligned_diel_inds(1);
for ii=num_cond+2:num_cond+num_diel
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+x_aligned_diel_inds(ii-num_cond);
end

if x_aligned_interface_inds(1)==0
    inds_glob(num_cond+num_diel+1,1)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+1,2)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+2,1)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+2,2)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+3,1)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+3,2)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+4,1)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+4,2)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+5,1)=inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+5,2)=inds_glob(num_cond+num_diel,2);
else
    inds_glob(num_cond+num_diel+1,1)=1+inds_glob(num_cond+num_diel,2);
    inds_glob(num_cond+num_diel+1,2)=inds_glob(num_cond+num_diel,2)+x_aligned_interface_inds(1);
    for ii=num_cond+num_diel+2:num_cond+num_diel+num_interface
        inds_glob(ii,1)=1+inds_glob(ii-1,2);
        inds_glob(ii,2)=inds_glob(ii-1,2)+x_aligned_interface_inds(ii-num_cond+num_diel);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%   y-aligned  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inds_glob(num_cond+num_diel+num_interface+1,1)=1+inds_glob(num_cond+num_diel+num_interface,2);
inds_glob(num_cond+num_diel+num_interface+1,2)=inds_glob(num_cond+num_diel+num_interface,2)+y_aligned_cond_inds(1);
for ii=num_cond+num_diel+num_interface+2:2*num_cond+num_diel+num_interface
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+y_aligned_cond_inds(ii-(num_cond+num_diel+num_interface));
end
inds_glob(2*num_cond+num_diel+num_interface+1,1)=1+inds_glob(2*num_cond+num_diel+num_interface,2);
inds_glob(2*num_cond+num_diel+num_interface+1,2)=inds_glob(2*num_cond+num_diel+num_interface,2)+y_aligned_diel_inds(1);
for ii=2*num_cond+num_diel+num_interface+2:2*(num_cond+num_diel)+num_interface
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+y_aligned_diel_inds(ii-(2*num_cond+num_diel+num_interface));
end

if y_aligned_interface_inds(1)==0
    inds_glob(2*(num_cond+num_diel)+num_interface+1,1)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+1,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)+num_interface+2,1)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+2,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)+num_interface+3,1)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+3,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)+num_interface+4,1)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+4,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)+num_interface+5,1)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+5,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2);
else
    inds_glob(2*(num_cond+num_diel)+num_interface+1,1)=1+inds_glob(2*(num_cond+num_diel)++num_interface,2);
    inds_glob(2*(num_cond+num_diel)++num_interface+1,2)=inds_glob(2*(num_cond+num_diel)++num_interface,2)+y_aligned_interface_inds(1);
    for ii=2*(num_cond+num_diel)+num_interface+2:2*(num_cond+num_diel+num_interface)
        inds_glob(ii,1)=1+inds_glob(ii-1,2);
        inds_glob(ii,2)=inds_glob(ii-1,2)+y_aligned_interface_inds(ii-(2*(num_cond+num_diel)+num_interface));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%   z-aligned  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inds_glob(2*(num_cond+num_diel+num_interface)+1,1)=1+inds_glob(2*(num_cond+num_diel+num_interface),2);
inds_glob(2*(num_cond+num_diel+num_interface)+1,2)=inds_glob(2*(num_cond+num_diel+num_interface),2)+z_aligned_cond_inds(1);
for ii=2*(num_cond+num_diel+num_interface)+2:3*num_cond+2*(num_diel+num_interface)
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+z_aligned_cond_inds(ii-(2*(num_cond+num_diel+num_interface)));
end
inds_glob(3*num_cond+2*(num_diel+num_interface)+1,1)=1+inds_glob(3*num_cond+2*(num_diel+num_interface),2);
inds_glob(3*num_cond+2*(num_diel+num_interface)+1,2)=inds_glob(3*num_cond+2*(num_diel+num_interface),2)+z_aligned_diel_inds(1);
for ii=3*num_cond+2*(num_diel+num_interface)+2:3*(num_cond+num_diel)+2*num_interface
    inds_glob(ii,1)=1+inds_glob(ii-1,2);
    inds_glob(ii,2)=inds_glob(ii-1,2)+z_aligned_diel_inds(ii-(3*num_cond+2*(num_diel+num_interface)));
end

if z_aligned_interface_inds(1)>0
    inds_glob(3*num_cond+3*num_diel+2*num_interface+1,1)=1+inds_glob(3*num_cond+3*num_diel+2*num_interface,2);
    inds_glob(3*num_cond+3*num_diel+2*num_interface+1,2)=inds_glob(3*num_cond+3*num_diel+2*num_interface,2)+z_aligned_interface_inds(1);
    for ii=3*num_cond+3*num_diel+2*num_interface+2:3*(num_cond+num_diel+num_interface)
        inds_glob(ii,1)=1+inds_glob(ii-1,2);
        inds_glob(ii,2)=inds_glob(ii-1,2)+z_aligned_interface_inds(ii-(3*num_cond+3*num_diel+2*num_interface));
    end
end
% inds_glob(86:89,1)=inds_glob(85,2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tend = toc(tini_pre);
sim_preproc=tend;
disp(['Total time for generate panels ::: ' ,num2str(sim_preproc)]);

s=whos('geom_bndry_panels'); 
memory_bndry_panels=s.bytes/1024^2;
memory_stage_preprocessing=memory_bndry_panels


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%          precond        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if fl_precond==1
%     [precond,ids,seg]=implement_preconditioner(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);ord=[];
% elseif fl_precond==2
% %     [precond,ids]=implement_preconditioner_ver2n(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);I=1;seg=1;
%     [precond,ids,ord]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);I=1;seg=1;
% else 
%     ids=1;precond=1;seg=1;ord=1;
% end
num_layer=[num_cond;num_diel;num_interface];
if fl_precond==2
    [precond,ids,ord,ind_layer]=precond_multi_layer(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,1,inds_glob,num_layer);seg=1;
else
    ids=1;precond=1;seg=1;ord=1;seg=1;
end
s=whos('precond'); 
memory_precond=s.bytes/1024^2;

memory_stage_preconditioner=memory_bndry_panels+memory_precond
% return
% % ------------------------------------------------------------------------
% %                  Obtain RHS Vector
% % ------------------------------------------------------------------------
num_unk=size(geom_bndry_panels,1);

if (num_diel > 0) % diel + conductor
    V_RHS=zeros(num_unk,num_cond);
    counter=num_cond+num_diel+num_interface;
    for ii=1:num_cond
        V_RHS(inds_glob(ii,1):inds_glob(ii,2),ii)=dx^2;
        V_RHS(inds_glob(ii+counter,1):inds_glob(ii+counter,2),ii)=dx^2;
        V_RHS(inds_glob(ii+2*counter,1):inds_glob(ii+2*counter,2),ii)=dx^2;
    end
        
        
    
% % the following will be fixed
%     dum_vect=zeros(num_x_aligned_cond+num_y_aligned_cond+num_z_aligned_cond,num_cond);
%     
%     for kk=1:num_cond
%         dum_vect(find(geom_bndry_panels(:,6) == kk),kk) = dx^2;
%     end
%     
%     V_RHS=zeros(num_unk,num_cond);
%     
%     temp=[inds_glob(1,1):inds_glob(1,2)];
%     temp=[temp [inds_glob(num_diel+2,1):inds_glob(num_diel+2,2)] [inds_glob(2*(num_diel+1)+1,1):inds_glob(2*(num_diel+1)+1,2)]];
%     for kk=1:num_cond
%         V_RHS(temp,kk)=dum_vect(:,kk);
%     end
%     
%     clear dum_vect temp

end

V_RHS_O=V_RHS;



% % ------------------------------------------------------------------------
% %                  Obtain the Circulant Tensors
% % ------------------------------------------------------------------------
tini = tic;
if fl_filling_or_retrieval==0
    size_up_lim=max([L,M,N]);
    if size_up_lim <=50
        load('fN_prestored_data_50.mat'); % according to the structure size to choose prestored Toeplitz
    elseif size_up_lim>50 && size_up_lim <=100
        load('fN_prestored_data_100.mat'); % according to the structure size to choose prestored Toeplitz
    elseif size_up_lim>100 && size_up_lim <=200
        load('fN_prestored_data_200.mat'); % according to the structure size to choose prestored Toeplitz
    elseif size_up_lim>200 && size_up_lim <=300
        load('fN_prestored_data_300.mat'); % according to the structure size to choose prestored Toeplitz
    elseif size_up_lim>300 && size_up_lim <=400
        load('fN_prestored_data_400.mat'); % according to the structure size to choose prestored Toeplitz
    else
        load('fN_prestored_data_500.mat'); % according to the structure size to choose prestored Toeplitz
    end
    [fN_ch_all]=retrival_circulant1(dx,L,M,N,Eps_inout,tole,fl_Tucker_decomp,fN_prestored_data);
else
    [fN_ch_all]=generate_circulant_tensor_charge(dx,L,M,N,Eps_inout,tole,fl_Tucker_decomp);
end
tend = toc(tini);
sim_circulant=tend;
disp(['Total time for circulant tensor ::: ' ,num2str(sim_circulant)]);

% % ------------------------------------------------------------------------
% %                  Iterative Solutions
% % ------------------------------------------------------------------------
% %  Iterative solution w/ FFT accelerated method   


num=[num_cond num_diel num_interface];

% rhoa_ch_dens=rand(num_unk,1);
% fMV_ch=matvect_mult_charge_pecdiel_multi_diel1(rhoa_ch_dens, fN_ch_all,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,L,M,N,num,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond,ind_layer);
% return

fMV_ch=@(rhoa_ch_dens)matvect_mult_charge_pecdiel_multi_diel1(rhoa_ch_dens, fN_ch_all,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,L,M,N,num,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond,ind_layer);

q_charge_vect=zeros(num_unk,num_cond);

disp(['Iterative solution started ... '])
tini = tic;
for kk=1:num_cond %num_cond
    if (kk==1)
        tini_out = tic;
    end
    [q_charge_vect(:,kk), flag, relres, iter, resvec] = pgmres(@(rhoa_ch_dens)fMV_ch(rhoa_ch_dens),...
        V_RHS(:,kk), inner_it, tol, outer_it);

    if (kk==1)
        sim_time = toc(tini_out);
        sim_iter = length(resvec);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      right precond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fl_precond==2
    %sss=1;%1/((dx^2)*((epsa(1)+epsb)/((epsa(1)-epsb)*2*eps0)));
    for kk=1:size(ids,1)
        q_charge_vect(ids{kk},:)=precond{kk}*q_charge_vect(ids{kk},:);
    end
    if isempty(ord)==0
        for kk=1:size(ord,1)
            %q_charge_vect(ord{kk},:)=sss*q_charge_vect(ord{kk},:);
            for ii=1:size(ord{kk})
            q_charge_vect(ord{kk}(ii),:)=(1/ss(ind_layer{kk}(ii)))*q_charge_vect(ord{kk}(ii),:);
            end
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tend = toc(tini);
sim_total_iter=tend;
disp(['Total time for iterative solution ::: ' ,num2str(tend)]);
disp(['Done... Iterative solution'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsa=[2.6 2.6];
% % load aaaa
% % q_charge_vect=q_charge_vect1;
% tot_num=num_cond+num_diel+num_interface;
% for ii=1:num_cond
%     if mod(ii,2)==1
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(1);
%     else
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(2);
%     end
% end
% for ii=tot_num+1:tot_num+num_cond
%     if mod(ii,2)==1
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(2);
%     else
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(1);
%     end
% end
% for ii=2*tot_num+1:2*tot_num+num_cond
%     if mod(ii,2)==1
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(1);
%     else
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(2);
%     end
% end
% for ii=39:39
%       q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)=q_charge_vect(inds_glob(ii,1):inds_glob(ii,2),:)*epsa(1);
% end

C_mat = V_RHS_O'*(q_charge_vect*2.6)
C_mat=real(C_mat)
%save('multi15_solution_new.mat','C_mat','geom_bndry_panels','q_charge_vect','dx','inds_glob')
%return

%%%---------------------------------------------------------------------
%                  Saving the data for post-processing
%%%---------------------------------------------------------------------
if (fl_save_data == 1)
    disp('-----------------------------------------------------')
    disp(['Saving Data...'])
    
%    data pertinent to geometry
    save('results_numexD_multilayer/data_geo.mat', 'dx', 'geom_bndry_panels', 'inds_glob','Eps_inout');
    
 %   data pertinent to solution
    save('results_numexD_multilayer/data_solution.mat', 'q_charge_vect','V_RHS','C_mat');
    
    disp(['Done... Saving data'])
    disp('-----------------------------------------------------')
end

% ------------------------------------------------------------------------
%                  Analytical Formulae
% ------------------------------------------------------------------------
% eps0=8.854187817e-12;
% radius=Dims(1,1)/2;
% radius_die=Dims(2,1)/2;% for dielectric sphere
% C_ana=4*pi*eps0*radius % for conductor sphere
% C_ana=4*pi*eps0*epsa*(radius*radius_die)/(radius*epsa+(radius_die-radius)) % for dielectric sphere
% 


fl_write_fastcap_file=0;
% % ------------------------------------------------------------------------
% %                  Writing FastCap Files
% % ------------------------------------------------------------------------

% preparing files for the fastcap
if (fl_write_fastcap_file == 1)
    
  id_int=[num_cond_x num_cond_z];
    
    for mm=1:size(geom_bndry_panels3,1)
        
        fid2=fopen(['intcon',num2str(id_int(1)),'by',num2str(id_int(2)),'_',num2str(mm),'.txt'],'w');
        
        fprintf(fid2,'%1s \n','0 Geometry file generated ...');
        fprintf(fid2,'%1s \n','* Cube with ... ');
        
        
        for kk=1:size(geom_bndry_panels3{mm},1)
            cen_pnl=geom_bndry_panels3{mm}(kk,1:3);
            if (geom_bndry_panels3{mm}(kk,4) == 1) % x-directed
                v1=[cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)-dx/2];
                v2=[cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)-dx/2];
                v3=[cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)+dx/2];
                v4=[cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)+dx/2];
            elseif(geom_bndry_panels3{mm}(kk,4) == -1) % -x-directed
                v1=[cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)+dx/2];
                v2=[cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)+dx/2];
                v3=[cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)-dx/2];
                v4=[cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)-dx/2];
            elseif(geom_bndry_panels3{mm}(kk,4)== 2) % y-directed
                v1=[cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)-dx/2];
                v2=[cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)-dx/2];
                v3=[cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)+dx/2];
                v4=[cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)+dx/2];
            elseif(geom_bndry_panels3{mm}(kk,4) == -2) % -y-directed
                v1=[cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)-dx/2];
                v2=[cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)-dx/2];
                v3=[cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)+dx/2];
                v4=[cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)+dx/2];
            elseif(geom_bndry_panels3{mm}(kk,4) == 3) % z-directed
                v1=[cen_pnl(1)+dx/2 cen_pnl(2)-dx/2 cen_pnl(3)];
                v2=[cen_pnl(1)+dx/2 cen_pnl(2)+dx/2 cen_pnl(3)];
                v3=[cen_pnl(1)-dx/2 cen_pnl(2)+dx/2 cen_pnl(3)];
                v4=[cen_pnl(1)-dx/2 cen_pnl(2)-dx/2 cen_pnl(3)];
            elseif(geom_bndry_panels3{mm}(kk,4) == -3) % -z-directed
                v1=[cen_pnl(1)-dx/2 cen_pnl(2)-dx/2 cen_pnl(3)];
                v2=[cen_pnl(1)-dx/2 cen_pnl(2)+dx/2 cen_pnl(3)];
                v3=[cen_pnl(1)+dx/2 cen_pnl(2)+dx/2 cen_pnl(3)];
                v4=[cen_pnl(1)+dx/2 cen_pnl(2)-dx/2 cen_pnl(3)];
            end
            
            fprintf(fid2,'%s %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n','Q 1 ',v1(1),v1(2),v1(3), v2(1),v2(2),v2(3), v3(1),v3(2),v3(3), v4(1),v4(2),v4(3) );
            
        end
        
        fclose(fid2);

    end
    disp('Done!...')
    
    % writing the main file for execution
    
    fid2=fopen(['intcon',num2str(id_int),'mm_dist.lst'],'w');
    fprintf(fid2,'%s\n','* By Gauss law, the capacitance is 4*pi*e0*er*a*b/(b-a(1-er)) ');
    fprintf(fid2,'%s\n','* Just some writing at the beginning, now actual lines are coming!!! ');
    
    for mm=1:size(geom_bndry_panels3,1)
        if (mm == size(geom_bndry_panels3,1))
            fl_name=['intcon',num2str(id_int(1)),'by',num2str(id_int(2)),'_',num2str(mm),'.txt'];
            all_string=['D ',fl_name,' 1.0 2.6  0.0 0.0 0.0 0.5 0.5 0.5 - '];
            fprintf(fid2,'%s\n',all_string);
            
        else
            fl_name=['intcon',num2str(id_int(1)),'by',num2str(id_int(2)),'_',num2str(mm),'.txt'];
            all_string=['C ',fl_name,' 2.6  0.0 0.0 0.0'];
            fprintf(fid2,'%s\n',all_string);
        end
    end
    
    
    fclose(fid2);
    
    disp('Done!...')  
    
    
end

% Execute postprocessor
post_processor_numexD