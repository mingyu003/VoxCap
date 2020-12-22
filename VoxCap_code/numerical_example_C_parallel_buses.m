clc; close all; clear all; format long e;

% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------
pre_define_the_path_for_folders
% -------------------------------------------------------------------------
%                  Inputs for Simulation
% -------------------------------------------------------------------------
epsa=[7];% dielectric permittivity, for conductor case, set it to null []
epsb=1;eps0=8.854187817e-12;
er = 0;  % epsilon_r of conductors
se=5.8e7; % conductivity of conductors
inner_it = 50; outer_it = 10; tol=1e-4; % iterative solver inputs
Res = 5 ; % voxel size (deltax)
tole=1e-8;%tolerance for Tucker decompression
fl_check_domain=0; % set to 1 for only plotting the structure
fl_check_geo=0; % set to 1 for only plotting the domain
fl_check_domain_pnls=0; % set to 1 for boundary panels of comp domain
fl_check_geo_pnls=0; % set to 1  for boundary panels of geometry
fl_save_data = 1; % set 1 for saving data for post-processing later on
fl_write_fastcap_file = 0; % set to 1 for outputing the panel info for fastcap sim
fl_Tucker_decomp=0; % set 1 to compress tensor with Tucker
fl_precond=2; %set 1 to implement preconditioner, 2 to preconditioner_ver2_n
num_vox_in_blk=10; % # of voxels in each preconditioner block. 
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
num_cond_x=5;
num_cond_z=12;

num_cond=num_cond_x*num_cond_z;
num_diel=length(epsa);

% cross-sections and lengths of buses
width_bus=80;
height_bus=20;
len_bus=800;
%len_bus=num_cond_x*width_bus
dist_y=width_bus+height_bus;
dist_z=2*height_bus; % between centers of 2 buses

if num_diel>0
    buff_btw_cond_diel=height_bus; % buffer to put dielectric around the conductors
    shft_vect=[buff_btw_cond_diel buff_btw_cond_diel buff_btw_cond_diel];
else
    shft_vect=[0 0 0];
end
% inputs for conductor
cen_cond=zeros(num_cond,3);
er_inout_cond=zeros(num_cond,3);
counter=1;
for kk=1:num_cond_z
    for ll=1:num_cond_x
%         if (mod(kk,2) == 1)
%             cen_cond(counter,1:3)=[width_bus/2+(dist*(ll-1)) len_bus/2 height_bus/2+(dist*(kk-1))] + shft_vect;
%         else
            cen_cond(counter,1:3)=[len_bus/2 width_bus/2+(dist_y*(ll-1)) height_bus/2+(dist_z*(kk-1))] + num_diel*shft_vect;
%         end
        if num_diel>0
            er_inout_cond(counter,1:2)=[er epsa(1)];
        else
            er_inout_cond(counter,1:2)=[er 0];
        end
        counter=counter+1;
    end
end

tmp = [len_bus width_bus height_bus];
counter=1;
Cnt_cond=[];Ref_pnts_cond=[];Eps_inout_cond=[];Dims_cond=[];Orients_cond=[];
for kk=1:num_cond_z
    for ll=1:num_cond_x
        Cnt_cond = [Cnt_cond; cen_cond(counter,1:3)];
        Ref_pnts_cond = [Ref_pnts_cond; cen_cond(counter,1:3)];
        Eps_inout_cond = [Eps_inout_cond; er_inout_cond(counter,1:2)];
        
        Dims_cond = [Dims_cond; tmp];
%         if (mod(kk,2) == 1)
%             Orients_cond=[Orients_cond;'y'];
%         else
            Orients_cond=[Orients_cond;'x'];
%         end
        counter=counter+1;
    end
end

% inputs for dielectric
if num_diel>0
    Dims_diel=zeros(num_diel,3);Cnt_diel=zeros(num_diel,3);Eps_inout_diel=zeros(num_diel,2);Orients_diel=[];
    for ii=1:num_diel
        Dims_diel(ii,1)=len_bus+2*ii*shft_vect(1); % length of dielectric_ii
        Dims_diel(ii,2)=width_bus+dist_y*(num_cond_x-1)+2*ii*shft_vect(2); % width of dielectric-ii
        Dims_diel(ii,3)=height_bus+dist_z*(num_cond_z-1)+2*ii*shft_vect(3); % height of dielectric-ii
        Cnt_diel(ii,:)=[Dims_diel(ii,1)/2  Dims_diel(ii,2)/2  Dims_diel(ii,3)/2]+ (num_diel-ii)*shft_vect;
        Orients_diel=[Orients_diel;'x'];
        if ii < num_diel
            Eps_inout_diel(ii,:)=[epsa(ii) epsa(ii+1)];
        else
            Eps_inout_diel(ii,:)=[epsa(ii) epsb];
        end
    end
end
if num_diel>0
    Ref_pnts_diel=Cnt_diel;
    Cnt=[Cnt_cond;Cnt_diel];
    Ref_pnts=[Ref_pnts_cond;Ref_pnts_diel];
    Eps_inout=[Eps_inout_cond;Eps_inout_diel];
else
    Cnt=Cnt_cond;
    Ref_pnts=Cnt_cond;
    Eps_inout=[Eps_inout_cond];
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
%bbox_max=[width_bus+dist*(num_cond_x-1) len_bus height_bus+dist*(num_cond_z-1)]+1e-8 ; % max coordinates of bbox
bbox_max=[len_bus width_bus+dist_y*(num_cond_x-1) height_bus+dist_z*(num_cond_z-1)]+1e-8 ; % max coordinates of bbox
if num_diel>0
    bbox_max=bbox_max+2*num_diel*shft_vect;
end
% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------
pre_print_out_inputs_generate_consts
timi1=tic;

% generate domain 3D grid
tic;
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);
disp(['Time for generating r ::: ',num2str(toc)])
[L,M,N,~] = size(r); % domain size
dx = Res; % voxel size
tic;
[idx,grid_intcon] = mingyu_distinguishvoxelsofstructure(dx,r,Cnt,Dims,Orients);
disp(['Time for generating grid_intcon ::: ',num2str(toc)])
clear Orients Dims;
% ------------------------------------------------------------------------
%                  Obtain Panel Coordinates and IDs
% -------------------------------------------------------------------------

disp('Generating panels of voxelized computational domain ...')
grid_tmp = ones(L,M,N); idx_dum = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp;
% idx_dum_cond_no=ones(length(idx_dum),1);
tic
[comp_dom_panels]=new_computational_panels(dx,r);
disp(['Time for generating computational domain ::: ',num2str(toc)])
%------------------------------------------------------------------------------
clear r idx_dum;% idx_dum_cond_no;
disp('Generating panels of voxelized individual bodies ...')

tic
geom_bndry_panels=cell(size(grid_intcon,1),1);

for kk=1:size(grid_intcon,1)
    idx_cd_no(1:length(idx{kk}))=kk;
    [geom_bndry_panels{kk}]=new_geom_panels(dx,grid_intcon{kk},idx{kk},idx_cd_no,num_vox_in_blk);
    clear idx_cd_no
end

clear grid_intcon idx;
disp(['Time for generating boundary panels ::: ',num2str(toc)])
s=whos('geom_bndry_panels'); 
memory_bndry_panels=s.bytes/1024^2;

if fl_check_geo_pnls>0
    if fl_check_geo_pnls==1
    plot_panel=cell2mat(geom_bndry_panels); 
%     else
%     plot_panel=cell2mat(geom_panels);
    end
    plot_panels(dx,plot_panel);
    clear plot_panel;
end
if fl_check_domain_pnls==1
    plot_panels(dx,comp_dom_panels);
end
%--------------------------------------------------------------------
% combining conductor panels
geom_bndry_panels3=geom_bndry_panels;

for ii=1:num_cond
    %     geom_bndry_panels{ii}(:,7)=0;
    geom_bndry_panels{ii}(:,8)=Eps_inout(1,1);
    geom_bndry_panels{ii}(:,9)=Eps_inout(1,2);
end
if num_diel>0
    for ii=1:num_diel
        geom_bndry_panels{num_cond+ii}(:,7)=ii;
        geom_bndry_panels{num_cond+ii}(:,8)=Eps_inout(num_cond+ii,1);
        geom_bndry_panels{num_cond+ii}(:,9)=Eps_inout(num_cond+ii,2);
    end
end


fl_plot_pnl_orients=0;
geom_bndry_panels=mingyu_pre_correct_pnl_orients(num_cond,num_diel,geom_bndry_panels,Eps_inout,Cnt,dx,fl_plot_pnl_orients);

if (num_diel == 0)
    fl_st_diel=1;
else
    fl_st_diel=num_diel+1;
end

num_panels=cell(fl_st_diel,1);
ids_panels=cell(fl_st_diel,1);
num_tot_panels=0;
geom_bndry_panels_cond=cell2mat(geom_bndry_panels(1:num_cond));
geom_bndry_panels_cond=sortrows(geom_bndry_panels_cond,5);
for kk=1:1
    [num_panels{kk},ids_panels{kk}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels_cond);
    num_tot_panels=num_tot_panels+num_panels{kk};
end

if num_diel>0
    for kk=1:num_diel
        [num_panels{kk+1},ids_panels{kk+1}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels{kk+num_cond});
        num_tot_panels=num_tot_panels+num_panels{kk+1};
    end
end
clear comp_dom_panels 

if (num_diel == 0)
    ids_panels{2}=cell(6,1);
end

geom_bndry_panels2=cell(fl_st_diel,1);
geom_bndry_panels2{1}=geom_bndry_panels_cond;
clear geom_bndry_panels_cond

if(num_diel > 0)% just add at the end
    for ii=1:num_diel
        geom_bndry_panels2{ii+1} = geom_bndry_panels{num_cond+ii};
    end
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
diel_pnl_id_locs=cell(num_diel,1);diel_pnl_ids=cell(num_diel,1);
for ii=1:num_diel
    diel_pnl_id_locs{ii}=cell(6,1); diel_pnl_ids{ii}=cell(6,1);
    if (num_diel > 0) % currently num_diel could be max 1
        % information for +/- direction directed panels from combined ids_panels
        diel_pnl_id_locs{ii}{1}=find(ids_panels{1+ii}{4}(:)>0);
        diel_pnl_ids{ii}{1}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{1});
        diel_pnl_id_locs{ii}{4}=find(ids_panels{1+ii}{4}(:)<0);
        diel_pnl_ids{ii}{4}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{4});
        
        diel_pnl_id_locs{ii}{2}=find(ids_panels{1+ii}{5}(:)>0);
        diel_pnl_ids{ii}{2}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{2});
        diel_pnl_id_locs{ii}{5}=find(ids_panels{1+ii}{5}(:)<0);
        diel_pnl_ids{ii}{5}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{5});
        
        diel_pnl_id_locs{ii}{3}=find(ids_panels{1+ii}{6}(:)>0);
        diel_pnl_ids{ii}{3}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{3});
        diel_pnl_id_locs{ii}{6}=find(ids_panels{1+ii}{6}(:)<0);
        diel_pnl_ids{ii}{6}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{6});
        
        % diel_pnl_id_locs stores the locations of elements of diel_pnl_ids which
        % are positive/negative x, y, and z directed panels.
        % Its entries {1} +x-directed, {2} +y-directed, {3} +z-directed, {4}
        % -x-directed, {5} -y directed, and {6} -z-directed
    end
end

% Attention: the following lines were added to distinguish x-, y-, z- aligned
% conductor and dielectric panels

num_x_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 1)); num_y_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 2)); num_z_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 3));
num_x_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0));
num_y_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0));
num_z_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0));


if num_diel>0
    num_x_aligned_diel=zeros(num_diel,1);num_y_aligned_diel=zeros(num_diel,1);num_z_aligned_diel=zeros(num_diel,1);
    for ii=1:num_diel
        num_x_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == ii));
        num_y_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == ii));
        num_z_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == ii));
    end
end

inds_glob=zeros(3*(num_diel+1),2);
inds_glob(1,:)=[1 num_x_aligned_cond];
inds_glob(num_diel+2,:) = [num_x_aligned_all+1 num_x_aligned_all+num_y_aligned_cond];
inds_glob(num_diel*2+3,:) = [num_x_aligned_all+num_y_aligned_all+1 num_x_aligned_all+num_y_aligned_all+num_z_aligned_cond]; % conductor part

if num_diel>0
    for ii=1:num_diel
        inds_glob(1+ii,:)=[inds_glob(ii,2)+1 inds_glob(ii,2)+sum(num_x_aligned_diel(ii))];
        inds_glob(num_diel+2+ii,:)=[inds_glob(num_diel+2+ii-1,2)+1 inds_glob(num_diel+2+ii-1,2)+sum(num_y_aligned_diel(ii))];
        inds_glob(num_diel*2+3+ii,:)=[inds_glob(num_diel*2+3+ii-1,2)+1 inds_glob(num_diel*2+3+ii-1,2)+sum(num_z_aligned_diel(ii))];
    end
end


pre_processing_time=toc(timi1)
% % ------------------------------------------------------------------------
% %                  Obtain RHS Vector
% % ------------------------------------------------------------------------
num_unk=size(geom_bndry_panels,1);
if (num_diel > 0) % diel + conductor
    % the following will be fixed
    dum_vect=zeros(num_x_aligned_cond+num_y_aligned_cond+num_z_aligned_cond,num_cond);
    
    for kk=1:num_cond
        dum_vect(find(geom_bndry_panels(:,6) == kk),kk) = dx^2;
    end
    
    V_RHS=zeros(num_unk,num_cond);
    
    temp=[inds_glob(1,1):inds_glob(1,2)];
    temp=[temp [inds_glob(num_diel+2,1):inds_glob(num_diel+2,2)] [inds_glob(2*(num_diel+1)+1,1):inds_glob(2*(num_diel+1)+1,2)]];
    for kk=1:num_cond
        V_RHS(temp,kk)=dum_vect(:,kk);
    end
    
    clear dum_vect temp
    
else % only conductor
    
    V_RHS=zeros(num_unk,num_cond);
    for kk=1:num_cond
        inds=find(abs(geom_bndry_panels(:,6)) == kk);
        V_RHS(inds,kk)=dx^2; % for only conductor panels
    end
    
    %V_RHS([[inds_glob(1,1):inds_glob(1,2)] [inds_glob(2,1):inds_glob(2,2)]  ...
    %    [inds_glob(3,1):inds_glob(3,2)]])=dx^2; % for only conductor panels
end
V_RHS_O=V_RHS;


if fl_precond==1
    [precond,ids,seg]=implement_preconditioner(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);
elseif fl_precond==2
%     [precond,ids]=implement_preconditioner_ver2n(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);seg=1;ord=[];
    [precond,ids,ord]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel);I=1;seg=1;
else 
    ids=1;precond=1;seg=1;ord=1;
end
% if fl_save_data==0
%     clear geom_bndry_panels
% end

s=whos('precond'); 
memory_precond=s.bytes/1024^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      left precond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if fl_precond==1
%     for ll=1:size(ids,1)
%         V_RHS(ids{ll},:)=precond{ll}*V_RHS(ids{ll},:);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------------------------------------------------------------------
% %                  Obtain the Circulant Tensors
% % ------------------------------------------------------------------------
tini=tic;
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
    [fN_ch_all]=retrival_circulant(dx,L,M,N,Eps_inout,tole,fl_Tucker_decomp,fN_prestored_data);
else
    [fN_ch_all]=generate_circulant_tensor_charge(dx,L,M,N,Eps_inout,tole,fl_Tucker_decomp);
end

if num_diel>0
    ss=zeros(num_diel,1);
    for ii=1:num_diel
        if ii < num_diel
            ss(ii)=(dx^2)*((epsa(ii)+epsa(ii+1))/((epsa(ii)-epsa(ii+1))*2*eps0)); % self-term for dielectric matrix
        else
            ss(ii)=(dx^2)*((epsa(ii)+epsb)/((epsa(ii)-epsb)*2*eps0)); % self-term for dielectric matrix
        end
    end
else
    ss=0;
end
tend = toc(tini);
sim_circulant=tend;
disp(['Total time for circulant tensor ::: ' ,num2str(sim_circulant)]);
% % ------------------------------------------------------------------------
% %                  Iterative Solutions
% % ------------------------------------------------------------------------

%  Iterative solution w/ FFT accelerated method
%--------------------------------------------------------------------
fl_multi_CPU=1;
if fl_multi_CPU==0
     warning off MATLAB:maxNumCompThreads:Deprecated
        maxNumCompThreads(1);
        LASTN = maxNumCompThreads;
end

% rhoa_ch_dens=rand(num_unk,1);
% fMV_ch=matvect_mult_charge_pecdiel_multi_diel(rhoa_ch_dens, fN_ch_all,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,L,M,N,num_diel,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond);
% 
% return

if num_diel>1
    fMV_ch=@(rhoa_ch_dens)matvect_mult_charge_pecdiel_multi_diel(rhoa_ch_dens, fN_ch_all,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,L,M,N,num_diel,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond);
else
    if num_diel==1
        diel_pnl_id_locs = diel_pnl_id_locs{1};  diel_pnl_ids = diel_pnl_ids{1};
    end
    fMV_ch=@(rhoa_ch_dens)matvect_mult_charge_pecdiel(rhoa_ch_dens, fN_ch_all,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,L,M,N,num_diel,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond);
end
%--------------------------------------------------------------------------------
q_charge_vect=zeros(num_unk,num_cond);

disp(['Iterative solution started ... '])
tini = tic;
for kk=1:num_cond
    [q_charge_vect(:,kk), flag, relres, iter, resvec] = pgmres(@(rhoa_ch_dens)fMV_ch(rhoa_ch_dens),...
        V_RHS(:,kk), inner_it, tol, outer_it);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      right precond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sss=1/((dx^2)*((epsa(1)+epsb)/((epsa(1)-epsb)*2*eps0)));
if fl_precond==2
    for kk=1:size(ids,1)
        q_charge_vect(ids{kk},:)=precond{kk}*q_charge_vect(ids{kk},:);
    end
    if isempty(ord)==0
        for kk=1:size(ord,1)
            q_charge_vect(ord{kk},:)=sss*q_charge_vect(ord{kk},:);
        end
    end
end

tend = toc(tini);
disp(['Total time for iterative solution ::: ' ,num2str(tend)]);
disp(['Done... Iterative solution'])

q_charge_vect=real(q_charge_vect);
% q_charge_vect=q_charge_vect(:,1);

% savefile = 'charge.mat';
% p = q_charge_vect;
% q = inds_glob;
% m=geom_bndry_panels;
% n=dx;
% save(savefile, 'p', 'q','m','n');
totalcharge_it_fast = sum(sum(q_charge_vect))

C_mat = V_RHS'*q_charge_vect;
C_mat=real(C_mat);
if (num_diel > 0)
    C_mat = V_RHS'*(q_charge_vect*(Eps_inout(1,2)));
end

% Saving the data for post-processing

if (fl_save_data == 1)
    disp('-----------------------------------------------------')
    disp(['Saving Data...'])
    
    % data pertinent to geometry
    %    data pertinent to geometry
    save('results_numexC_parallel_buses/data_geo.mat', 'dx', 'geom_bndry_panels', 'inds_glob','Eps_inout');
    
 %   data pertinent to solution
    save('results_numexC_parallel_buses/data_solution.mat', 'q_charge_vect','V_RHS','C_mat');
    
    disp(['Done... Saving data'])
    disp('-----------------------------------------------------')
end


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
            all_string=['D ',fl_name,' 1.0 7.0  0.0 0.0 0.0 220 260 250 - '];
            fprintf(fid2,'%s\n',all_string);
            
        else
            fl_name=['intcon',num2str(id_int(1)),'by',num2str(id_int(2)),'_',num2str(mm),'.txt'];
            all_string=['C ',fl_name,' 7.0  0.0 0.0 0.0'];
            fprintf(fid2,'%s\n',all_string);
        end
    end
    
    
    fclose(fid2);
    
    disp('Done!...')
    
end
% Execute postprocessor
post_processor_numexC
