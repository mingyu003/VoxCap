function [z_mat] = slow_method_z_mat(dx,geom_bndry_panels,fl_full_or_blocks)
% this function is used to get system matrix
%%%INPUTS
%dx: resolution, the border length of each voxel
%fl_full_or_blocks:0: partition Z_mat into blocks, 1: use full Z_mat
%geom_bndry_panels: there are 9 columns.
%   columns 1-3: center point coordinates of boundary panels.
%   column 4: the normal direction of panels
%   column 5: the numbering of panels
%   column 6: which conductor or dielectric the panel belong to.
%   column 7: flag, 0 represents conductor panel, 1 represents dielectrical panel
%   columns 8-9: inner and outer permittivity
%%%OTPUT
%z_mat: system matrix

disp('-----------------------------------------------------')
if (fl_full_or_blocks == 1)
    disp('Computing full system matrix')
else
    disp('Computing full system matrix blocks')
end

 fl_diel=0;
if (sum(geom_bndry_panels(:,7)) > 1e-13); % check whether we have a dielectric panel
    fl_diel = 1; % yes, we have
end

%-------------------------------------------------------------------------
%                  numenrical quadrature parameters preparation
%-------------------------------------------------------------------------
% Quadrature will be used in the hybrid calling
num_quad_pnts = 3; % set to 7 if wanted unnecessarily more accuracy
%num_quad_pnts: Gauss-quadrature sample points numbers
num_diff_pnts = 2; % set to 4 if wanted unnecessarily more accuracy
%num_diff_pnts: when using numenrical partial differentiation, use it to
%get finite differentiation coefficients
[smpl_pnts,smpl_wghts] = weights_points(num_quad_pnts,4,[-1 1;-1 1;-1 1;-1 1]);
smpl_wghts = smpl_wghts*(dx^4)/16; smpl_pnts = smpl_pnts*(dx)/2;
%--------------------------------------------------------------------------
if (fl_full_or_blocks == 1)
    num_unk=size(geom_bndry_panels,1);
    disp(['# of unknowns for system matrix :::',num2str(num_unk)])
    z_mat=zeros(num_unk,num_unk);
    tic
    parfor kk=1:num_unk % observer panel
        %if (mod(kk,20) == 0); disp([num2str(num_unk),'  ',num2str(kk)]); end;
        for ll=1:num_unk % source panel
%             if geom_bndry_panels(kk,7)>0 && kk==ll
%                 z_mat(kk,ll)=0;
%             else
            [z_mat(kk,ll)]=compute_Galerkin_int_cond_diel_integral(dx,geom_bndry_panels(ll,1:3),geom_bndry_panels(ll,4), ...
                geom_bndry_panels(kk,1:3),geom_bndry_panels(kk,4),geom_bndry_panels(kk,7), ...
                geom_bndry_panels(kk,8),geom_bndry_panels(kk,9),smpl_pnts,smpl_wghts,num_diff_pnts);
%             end
        end
    end
    disp(['Time for computing system matrix ::: ',num2str(toc)])
    eps0=8.854187817e-12;
    one_over_4pieps0 = (1/(4*pi*eps0));
    z_mat = z_mat * one_over_4pieps0;
    
    if (fl_diel == 1)
        
        % The format of Z_mat for conductor+dielectric case is for multiplication
        % with the charge coefficient vector in the following format
        % [rhoa_x_c rhoa_y_c rhoa_z_c rhoa_x_d rhoa_y_d rhoa_z_d]
        
        % To have the Z_mat for multiplication with the charge coefficient vector
        % of format [rhoa_x_c rhoa_x_d rhoa_y_c rhoa_y_d rhoa_z_c rhoa_z_d],
        % the following additions are included:
        
        % Lets find out x-, y-, anz z- aligned dielectric/conductor/all panels
        x_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0);
        x_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 1);
        
        y_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0);
        y_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 1);
        
        z_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0);
        z_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 1);
        
        inds_sorted=[x_aligned_cond_inds; x_aligned_diel_inds; ...
        y_aligned_cond_inds; y_aligned_diel_inds; z_aligned_cond_inds; z_aligned_diel_inds];
%         x_aligned_cond_inds1=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 1);
% x_aligned_cond_inds2=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 2);
% x_aligned_cond_inds3=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 3);
% x_aligned_cond_inds4=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 4);
% x_aligned_cond_inds5=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 5);
% x_aligned_cond_inds6=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 6);
% x_aligned_cond_inds7=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 7);
% x_aligned_cond_inds8=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 8);
% x_aligned_cond_inds9=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 9);
% x_aligned_cond_inds10=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 10);
% x_aligned_diel_inds1=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 11);
% x_aligned_diel_inds2=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 12);
% 
% y_aligned_cond_inds1=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 1);
% y_aligned_cond_inds2=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 2);
% y_aligned_cond_inds3=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 3);
% y_aligned_cond_inds4=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 4);
% y_aligned_cond_inds5=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 5);
% y_aligned_cond_inds6=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 6);
% y_aligned_cond_inds7=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 7);
% y_aligned_cond_inds8=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 8);
% y_aligned_cond_inds9=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 9);
% y_aligned_cond_inds10=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 10);
% y_aligned_diel_inds1=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 11);
% y_aligned_diel_inds2=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 12);
% 
% z_aligned_cond_inds1=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 1);
% z_aligned_cond_inds2=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 2);
% z_aligned_cond_inds3=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 3);
% z_aligned_cond_inds4=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 4);
% z_aligned_cond_inds5=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 5);
% z_aligned_cond_inds6=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 6);
% z_aligned_cond_inds7=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 7);
% z_aligned_cond_inds8=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 8);
% z_aligned_cond_inds9=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 9);
% z_aligned_cond_inds10=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == 10);
% z_aligned_diel_inds1=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 11);
% z_aligned_diel_inds2=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) >= 1 & geom_bndry_panels(:,6) == 12);
% 
% 
% inds_sorted=[x_aligned_cond_inds1;x_aligned_cond_inds2;x_aligned_cond_inds3;x_aligned_cond_inds4;x_aligned_cond_inds5;x_aligned_cond_inds6; ...
%    x_aligned_cond_inds7;x_aligned_cond_inds8;x_aligned_cond_inds9;x_aligned_cond_inds10;x_aligned_diel_inds1;x_aligned_diel_inds2; ...
%     y_aligned_cond_inds1;y_aligned_cond_inds2;y_aligned_cond_inds3;y_aligned_cond_inds4;y_aligned_cond_inds5;y_aligned_cond_inds6; ...
%     y_aligned_cond_inds7;y_aligned_cond_inds8;y_aligned_cond_inds9;y_aligned_cond_inds10; y_aligned_diel_inds1;y_aligned_diel_inds2; ...
%     z_aligned_cond_inds1;z_aligned_cond_inds2;z_aligned_cond_inds3;z_aligned_cond_inds4;z_aligned_cond_inds5;z_aligned_cond_inds6; ...
%     z_aligned_cond_inds7;z_aligned_cond_inds8;z_aligned_cond_inds9;z_aligned_cond_inds10;z_aligned_diel_inds1; z_aligned_diel_inds2];


        z_mat=z_mat(inds_sorted,inds_sorted);
        
    end
else
    tic
    z_mat=cell(3,3,2);
    %figure out the which panels are conductor panels & dielectric panels
    %observer_conductor
    x_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 1);
    y_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 1);
    z_aligned_diel_inds=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 1);
    %observer_conductor
    x_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0);
    y_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0);
    z_aligned_cond_inds=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0);
    %source
    x_aligned_all_inds=find(abs(geom_bndry_panels(:,4)) == 1);
    y_aligned_all_inds=find(abs(geom_bndry_panels(:,4)) == 2);
    z_aligned_all_inds=find(abs(geom_bndry_panels(:,4)) == 3);
    
    % Obtain blocks
    % (:,:,1), for conductor observer; (:,:,2), for dielectric observer
    
    % Gxx for conductor and dielectric observer
    [z_mat{1,1,1}]=obtain_z_mat_block(dx,x_aligned_cond_inds,x_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{1,1,2}]=obtain_z_mat_block(dx,x_aligned_diel_inds,x_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    
    % Gyy for onductor and dielectric observer
    [z_mat{2,2,1}]=obtain_z_mat_block(dx,y_aligned_cond_inds,y_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{2,2,2}]=obtain_z_mat_block(dx,y_aligned_diel_inds,y_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    
    % Gzz for onductor and dielectric observer
    [z_mat{3,3,1}]=obtain_z_mat_block(dx,z_aligned_cond_inds,z_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{3,3,2}]=obtain_z_mat_block(dx,z_aligned_diel_inds,z_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    %%% Upper triangle
    % Gxy - - source y-aligned panel, observer x-aligned panel
    [z_mat{1,2,1}]=obtain_z_mat_block(dx,x_aligned_cond_inds,y_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{1,2,2}]=obtain_z_mat_block(dx,x_aligned_diel_inds,y_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    % Gxz - - source z-aligned panel, observer x-aligned panel
    [z_mat{1,3,1}]=obtain_z_mat_block(dx,x_aligned_cond_inds,z_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{1,3,2}]=obtain_z_mat_block(dx,x_aligned_diel_inds,z_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    % Gyz - - source z-aligned panel, observer y-aligned panel
    [z_mat{2,3,1}]=obtain_z_mat_block(dx,y_aligned_cond_inds,z_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{2,3,2}]=obtain_z_mat_block(dx,y_aligned_diel_inds,z_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    %%% Lower triangle
    % Gyx - - source x-aligned panel, observer y-aligned panel
    [z_mat{2,1,1}]=obtain_z_mat_block(dx,y_aligned_cond_inds,x_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{2,1,2}]=obtain_z_mat_block(dx,y_aligned_diel_inds,x_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    % Gzx - - source x-aligned panel, observer z-aligned panel
    [z_mat{3,1,1}]=obtain_z_mat_block(dx,z_aligned_cond_inds,x_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    
    if (fl_diel == 1)
        [z_mat{3,1,2}]=obtain_z_mat_block(dx,z_aligned_diel_inds,x_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    % Gzy - - source y-aligned panel, observer z-aligned panel
    [z_mat{3,2,1}]=obtain_z_mat_block(dx,z_aligned_cond_inds,y_aligned_all_inds,...
        geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    if (fl_diel == 1)
        [z_mat{3,2,2}]=obtain_z_mat_block(dx,z_aligned_diel_inds,y_aligned_all_inds,...
            geom_bndry_panels,smpl_pnts,smpl_wghts,num_diff_pnts);
    end
    if (fl_diel == 0)
        z_mat=cell2mat(z_mat);
    elseif (fl_diel == 1)
        z_mat={z_mat{1,1,1},z_mat{1,2,1},z_mat{1,3,1}; ...
               z_mat{1,1,2},z_mat{1,2,2},z_mat{1,3,2}; ...
               z_mat{2,1,1},z_mat{2,2,1},z_mat{2,3,1}; ...
               z_mat{2,1,2},z_mat{2,2,2},z_mat{2,3,2}; ...
               z_mat{3,1,1},z_mat{3,2,1},z_mat{3,3,1}; ...
               z_mat{3,1,2},z_mat{3,2,2},z_mat{3,3,2}};
        z_mat=cell2mat(z_mat);
    end
    
    disp(['Time for computing all system matrix blocks ::: ',num2str(toc)])
    
end

if (fl_full_or_blocks == 1)
    disp('Done... Computing system matrix ')
else
    disp('Done... Computing system matrix blocks ')
end
disp('-----------------------------------------------------')