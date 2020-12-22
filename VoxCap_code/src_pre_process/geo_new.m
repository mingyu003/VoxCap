function [bndry_panel,boolean_tens_cell]=geo_new(dx,x_bnd,y_bnd,z_bnd,grid_intcon,idxS,idxS_cond,fl_precond)
%%%%the function can just be used in convex structure!!!!!!!!!!!!!!!

[L,M,N,~]=size(grid_intcon);
%%%------------------------------------------------------------------
%             obtain geom panels' coordinates and normal
%%%------------------------------------------------------------------

% total_panels=zeros(L,M,N,24);
boolean_tens=zeros(L,M,N);
boolean_tens(idxS)=1;

ids_boolean=zeros(L*M*N,1);
for nn=2:N-1
    for mm=2:M-1
        for ll=2:L-1
            if boolean_tens(ll,mm,nn)==1 && boolean_tens(ll+1,mm,nn)==1 && boolean_tens(ll-1,mm,nn)==1 && boolean_tens(ll,mm+1,nn)==1 ...
                    && boolean_tens(ll,mm-1,nn)==1 && boolean_tens(ll,mm,nn+1)==1 && boolean_tens(ll,mm,nn-1)==1
                ids_boolean(ll+(mm-1)*L+(nn-1)*L*M)=1;
            end
        end
    end
end
boolean_tens(find(ids_boolean==1))=0;
if fl_precond==1
    boolean_tens_cell=cell(5,1);
    [boolean_tens_cell{1},boolean_tens_cell{2},boolean_tens_cell{3},boolean_tens_cell{4},boolean_tens_cell{5}] = Tucker(boolean_tens,1e-4);
else
    boolean_tens_cell=1;
end
% if fl_precond==1
%     boolean_tens_cell=cell(N,1);
%     for ii=1:N
%         boolean_tens_cell{ii}=sparse(boolean_tens(:,:,ii));
%     end
% else
%     boolean_tens_cell=1;
% end

dumxr=0; % x right
dumyb=0; % y back
dumzu=0; % z up
dumxl=0; % x left
dumyf=0; % y front
dumzb=0; % z bottom

bndry_x_r=zeros(L,12);
bndry_x_l=zeros(L,12);
bndry_y_f=zeros(M,12);
bndry_y_b=zeros(M,12);
bndry_z_b=zeros(N,12);
bndry_z_u=zeros(N,12);

% left, right panels, direction=1
for mm=1:M
    for nn=1:N
        ids_x=find(boolean_tens(:,mm,nn)==1);
        if length(ids_x)>1e-13
            ids_right=max(ids_x);
            ids_left=min(ids_x);
            dumxr=dumxr+1;
            dumxl=dumxl+1;
            bndry_x_r(dumxr,1)=grid_intcon(ids_right,mm,nn,1)+dx/2;
            bndry_x_r(dumxr,2)=grid_intcon(ids_right,mm,nn,2);
            bndry_x_r(dumxr,3)=grid_intcon(ids_right,mm,nn,3);
            bndry_x_r(dumxr,4)=1;
            bndry_x_r(dumxr,10:12)=[ids_right,mm,nn];
            bndry_x_l(dumxl,1)=grid_intcon(ids_left,mm,nn,1)-dx/2;
            bndry_x_l(dumxl,2)=grid_intcon(ids_left,mm,nn,2);
            bndry_x_l(dumxl,3)=grid_intcon(ids_left,mm,nn,3);
            bndry_x_l(dumxl,4)=1;
            bndry_x_l(dumxl,10:12)=[ids_left,mm,nn];
            clear ids_x
        end
    end
end

% front, back panels, direction=2
for ll=1:L
    for nn=1:N
        ids_y=find(boolean_tens(ll,:,nn)==1);
        if length(ids_y)>1e-13
            ids_back=max(ids_y);
            ids_front=min(ids_y);
            dumyb=dumyb+1;
            dumyf=dumyf+1;
            bndry_y_b(dumyb,1)=grid_intcon(ll,ids_back,nn,1);
            bndry_y_b(dumyb,2)=grid_intcon(ll,ids_back,nn,2)+dx/2;
            bndry_y_b(dumyb,3)=grid_intcon(ll,ids_back,nn,3);
            bndry_y_b(dumyb,4)=2;
            bndry_y_b(dumyb,10:12)=[ll,ids_back,nn];
            bndry_y_f(dumyf,1)=grid_intcon(ll,ids_front,nn,1);
            bndry_y_f(dumyf,2)=grid_intcon(ll,ids_front,nn,2)-dx/2;
            bndry_y_f(dumyf,3)=grid_intcon(ll,ids_front,nn,3);
            bndry_y_f(dumyf,4)=2;
            bndry_y_f(dumyf,10:12)=[ll,ids_front,nn];
            clear ids_y
        end
    end
end

% up, bottom panels, direction=3
for ll=1:L
    for mm=1:M
        ids_z=find(boolean_tens(ll,mm,:)==1);
        if length(ids_z)>1e-13
            ids_up=max(ids_z);
            ids_bottom=min(ids_z);
            dumzu=dumzu+1;
            dumzb=dumzb+1;
            bndry_z_u(dumzu,1)=grid_intcon(ll,mm,ids_up,1);
            bndry_z_u(dumzu,2)=grid_intcon(ll,mm,ids_up,2);
            bndry_z_u(dumzu,3)=grid_intcon(ll,mm,ids_up,3)+dx/2;
            bndry_z_u(dumzu,4)=3;
            bndry_z_u(dumzu,10:12)=[ll,mm,ids_up];
            bndry_z_b(dumzb,1)=grid_intcon(ll,mm,ids_bottom,1);
            bndry_z_b(dumzb,2)=grid_intcon(ll,mm,ids_bottom,2);
            bndry_z_b(dumzb,3)=grid_intcon(ll,mm,ids_bottom,3)-dx/2;
            bndry_z_b(dumzb,4)=3;
            bndry_z_b(dumzb,10:12)=[ll,mm,ids_bottom];
            clear ids_z
        end
    end
end

bndry_panel=[bndry_x_l; bndry_x_r; bndry_y_f; bndry_y_b; bndry_z_b; bndry_z_u];
bndry_panel(all(bndry_panel==0,2),:)=[];
%%%-------------------------------------------------------------------
%                  obtain numbering for geom panels
%%%-------------------------------------------------------------------
% generate ids of x_directed panels
tensor_ijk_xdir_pnl=zeros(L+1,M,N);
dum_cnt=1;
for mm=1:N% z-variation
    for ll=1:M % y-variation
        for kk=1:L+1 % x-variation
            tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end
% generate ids of  y_directed panels
tensor_ijk_ydir_pnl=zeros(L,M+1,N);
for mm=1:N % z-variation
    for ll=1:M+1 % y-variation
        for kk=1:L % x-variation
            tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end

% generate ids of z_directed panels
tensor_ijk_zdir_pnl=zeros(L,M,N+1);
for mm=1:N+1 % z-variation
    for ll=1:M % y-variation
        for kk=1:L % x-variation
            tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
            dum_cnt=dum_cnt+1;
        end
    end
end
% % generate ids of  y_directed panels
% tensor_ijk_ydir_pnl=zeros(L,M+1,N);
% for kk=1:L % x-variation
%     for mm=1:N % z-variation
%         for ll=1:M+1 % y-variation
%             tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
%             dum_cnt=dum_cnt+1;
%         end
%     end
% end
% 
% % generate ids of z_directed panels
% tensor_ijk_zdir_pnl=zeros(L,M,N+1);
% for ll=1:M % y-variation
%     for kk=1:L % x-variation
%         for mm=1:N+1 % z-variation
%             tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
%             dum_cnt=dum_cnt+1;
%         end
%     end
% end

bbox_origin=[min(x_bnd) min(y_bnd) min(z_bnd)];
org_x_panel=[bbox_origin(1) bbox_origin(2)+dx/2 bbox_origin(3)+dx/2];
org_y_panel=[bbox_origin(1)+dx/2 bbox_origin(2) bbox_origin(3)+dx/2];
org_z_panel=[bbox_origin(1)+dx/2 bbox_origin(2)+dx/2 bbox_origin(3)];

%%%---------------------------------------------------------------
%                    find bundry panels numbering
%%%---------------------------------------------------------------
for kk=1:size(bndry_panel,1)
    if (bndry_panel(kk,4) == 1) % x-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_x_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 2) % y-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_y_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 3) % z-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_z_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end
ids=idxS_cond(1);
bndry_panel(:,6)=ids;


disp('Done... Extracting unique panels of a grid')
disp('-----------------------------------------------------')