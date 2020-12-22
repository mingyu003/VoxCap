function [bndry_panel]=new_geom_panels(dx,grid_intcon,idxS,idxS_cond,num_vox_in_blk)

% the function is used to generate geometry panels and boundary panels
% fl_check_pnls=1, plot all geometry panels
% fl_check_pnls=2, plot boundary panels


[L,M,N,~]=size(grid_intcon);
num_blk_x=ceil(L/num_vox_in_blk);  num_blk_y=ceil(M/num_vox_in_blk);  num_blk_z=ceil(N/num_vox_in_blk);
bound=floor([L,M,N]./num_vox_in_blk)*num_vox_in_blk;

if mod(L,num_vox_in_blk)>0 && L>1
    num_big_blk_x=num_blk_x-1;
else
    num_big_blk_x=num_blk_x;
end
if mod(M,num_vox_in_blk)>0 && M>1
    num_big_blk_y=num_blk_y-1;
else
    num_big_blk_y=num_blk_y;
end
if mod(N,num_vox_in_blk)>0 && N>1
    num_big_blk_z=num_blk_z-1;
else
    num_big_blk_z=num_blk_z;
end

num_blk_x=num_big_blk_x;  num_blk_y=num_big_blk_y;  num_blk_z=num_big_blk_z;

%%%------------------------------------------------------------------
%             obtain geom panels' coordinates and normal
%%%------------------------------------------------------------------
boolean_tens=zeros(L,M,N);
boolean_tens(idxS)=1;
dumxr=0; % x right
dumyb=0; % y back
dumzu=0; % z up
dumxl=0; % x left
dumyf=0; % y front
dumzb=0; % z bottom

bndry_x_r=zeros(L*M*N,10);
bndry_x_l=zeros(L*M*N,10);
bndry_y_f=zeros(L*M*N,10);
bndry_y_b=zeros(L*M*N,10);
bndry_z_b=zeros(L*M*N,10);
bndry_z_u=zeros(L*M*N,10);
total_panels= total_panel(boolean_tens, grid_intcon,dx);
temp=zeros(1,3);

for kk=1:L
    for ll=1:M
        for mm=1:N
            if (boolean_tens(kk,ll,mm))==1%-1)<1e-12
                
                if ((kk+1)<=L && boolean_tens(kk+1,ll,mm)<1e-12) || kk==L  % right bndry
                    dumxr=dumxr+1;
                    bndry_x_r(dumxr,1:3)=total_panels(kk,ll,mm,5:7);
                    bndry_x_r(dumxr,4)=1;
%                     bndry_x_r(dumxr,10)=kk;
%                     bndry_x_r(dumxr,11)=ll;
%                     bndry_x_r(dumxr,12)=mm;
%                     bndry_x_r(dumxr,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_x_r(dumxr,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                
                if ((ll+1)<=M &&  boolean_tens(kk,ll+1,mm)<1e-12) || ll==M  %back bndry
                    
                    dumyb=dumyb+1;
                    bndry_y_b(dumyb,1:3)=total_panels(kk,ll,mm,13:15);
                    bndry_y_b(dumyb,4)=2;
%                     bndry_y_b(dumyb,10)=kk;
%                     bndry_y_b(dumyb,11)=ll;
%                     bndry_y_b(dumyb,12)=mm;
%                     bndry_y_b(dumyb,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_y_b(dumyb,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                if ((mm+1)<=N && boolean_tens(kk,ll,mm+1)<1e-12) || mm==N%up bndry
                    dumzu=dumzu+1;
                    bndry_z_u(dumzu,1:3)=total_panels(kk,ll,mm,21:23);
                    bndry_z_u(dumzu,4)=3;
%                     bndry_z_u(dumzu,10)=kk;
%                     bndry_z_u(dumzu,11)=ll;
%                     bndry_z_u(dumzu,12)=mm;
%                     bndry_z_u(dumzu,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_z_u(dumzu,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                if ((kk-1)>=1 && boolean_tens(kk-1,ll,mm)<1e-12) || kk==1 %left bndry
                    dumxl=dumxl+1;
                    bndry_x_l(dumxl,1:3)=total_panels(kk,ll,mm,1:3);
                    bndry_x_l(dumxl,4)=1;
%                     bndry_x_l(dumxl,10)=kk;
%                     bndry_x_l(dumxl,11)=ll;
%                     bndry_x_l(dumxl,12)=mm;
%                     bndry_x_l(dumxl,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_x_l(dumxl,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                if ((ll-1)>=1 && boolean_tens(kk,ll-1,mm)<1e-12)|| ll==1%front bndry
                    dumyf=dumyf+1;
                    bndry_y_f(dumyf,1:3)=total_panels(kk,ll,mm,9:11);
                    bndry_y_f(dumyf,4)=2;
%                     bndry_y_f(dumyf,10)=kk;
%                     bndry_y_f(dumyf,11)=ll;
%                     bndry_y_f(dumyf,12)=mm;
%                     bndry_y_f(dumyf,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_y_f(dumyf,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                if ((mm-1)>=1 && boolean_tens(kk,ll,mm-1)<1e-12) || mm==1 %bottom bndry
                    dumzb=dumzb+1;
                    bndry_z_b(dumzb,1:3)=total_panels(kk,ll,mm,17:19);
                    bndry_z_b(dumzb,4)=3;
%                     bndry_z_b(dumzb,10)=kk;
%                     bndry_z_b(dumzb,11)=ll;
%                     bndry_z_b(dumzb,12)=mm;
%                     bndry_z_b(dumzb,10)=kk+(ll-1)*L+(mm-1)*L*M;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
%                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_z_b(dumzb,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
            end
        end  
    end
end
% 
% for kk=1:L
%     for ll=1:M
%         for mm=1:N
%             if (boolean_tens(kk,ll,mm))==1%-1)<1e-12
% %                 total_panels(kk,ll,mm,1:3)=grid_intcon(kk,ll,mm,1:3);
% %                 total_panels(kk,ll,mm,9:11)=grid_intcon(kk,ll,mm,1:3);
% %                 total_panels(kk,ll,mm,17:19)=grid_intcon(kk,ll,mm,1:3);
% %                 total_panels(kk,ll,mm,1)=total_panels(kk,ll,mm,1)-dx/2;
% %                 total_panels(kk,ll,mm,10)=total_panels(kk,ll,mm,10)-dx/2;
% %                 total_panels(kk,ll,mm,19)=total_panels(kk,ll,mm,19)-dx/2;
% %                 total_panels(kk,ll,mm,4)=1;
% %                 total_panels(kk,ll,mm,12)=2;
% %                 total_panels(kk,ll,mm,20)=3;
%                 total_panels(kk,ll,mm,1)=grid_intcon(kk,ll,mm,1)-dx/2;
%                 total_panels(kk,ll,mm,2)=grid_intcon(kk,ll,mm,2);
%                 total_panels(kk,ll,mm,3)=grid_intcon(kk,ll,mm,3);
%                 total_panels(kk,ll,mm,4)=1;
%                 total_panels(kk,ll,mm,9)=grid_intcon(kk,ll,mm,1);
%                 total_panels(kk,ll,mm,10)=grid_intcon(kk,ll,mm,2)-dx/2;
%                 total_panels(kk,ll,mm,11)=grid_intcon(kk,ll,mm,3);
%                 total_panels(kk,ll,mm,12)=2;
%                 total_panels(kk,ll,mm,17)=grid_intcon(kk,ll,mm,1);
%                 total_panels(kk,ll,mm,18)=grid_intcon(kk,ll,mm,2);
%                 total_panels(kk,ll,mm,19)=grid_intcon(kk,ll,mm,3)-dx/2;
%                 total_panels(kk,ll,mm,20)=3;
%                 
%                 if ((kk+1)<=L && boolean_tens(kk+1,ll,mm)<1e-12) || kk==L  % right bndry
%                     dumxr=dumxr+1;
% %                     total_panels(kk,ll,mm,5:7)=grid_intcon(kk,ll,mm,1:3);
% %                     total_panels(kk,ll,mm,5)=total_panels(kk,ll,mm,5)+dx/2;
%                     total_panels(kk,ll,mm,5)=grid_intcon(kk,ll,mm,1)+dx/2;
%                     total_panels(kk,ll,mm,6)=grid_intcon(kk,ll,mm,2);
%                     total_panels(kk,ll,mm,7)=grid_intcon(kk,ll,mm,3);
%                     total_panels(kk,ll,mm,8)=1;
%                     bndry_x_r(dumxr,1:3)=total_panels(kk,ll,mm,5:7);
%                     bndry_x_r(dumxr,4)=1;
%                 end
%                 
%                 
%                 if ((ll+1)<=M &&  boolean_tens(kk,ll+1,mm)<1e-12) || ll==M  %back bndry
%                     
%                     dumyb=dumyb+1;
% %                     total_panels(kk,ll,mm,13:15)=grid_intcon(kk,ll,mm,1:3);
% %                     total_panels(kk,ll,mm,14)=total_panels(kk,ll,mm,14)+dx/2;
%                     total_panels(kk,ll,mm,13)=grid_intcon(kk,ll,mm,1);
%                     total_panels(kk,ll,mm,14)=grid_intcon(kk,ll,mm,2)+dx/2;
%                     total_panels(kk,ll,mm,15)=grid_intcon(kk,ll,mm,3);
%                     total_panels(kk,ll,mm,16)=2;
%                     bndry_y_b(dumyb,1:3)=total_panels(kk,ll,mm,13:15); 
%                     bndry_y_b(dumyb,4)=2;
%                 end
%                                
%                 if ((mm+1)<=N && boolean_tens(kk,ll,mm+1)<1e-12) || mm==N%up bndry                    
%                     dumzu=dumzu+1;
% %                     total_panels(kk,ll,mm,21:23)=grid_intcon(kk,ll,mm,1:3);
% %                     total_panels(kk,ll,mm,23)=total_panels(kk,ll,mm,23)+dx/2;
%                     total_panels(kk,ll,mm,21)=grid_intcon(kk,ll,mm,1);
%                     total_panels(kk,ll,mm,22)=grid_intcon(kk,ll,mm,2);
%                     total_panels(kk,ll,mm,23)=grid_intcon(kk,ll,mm,3)+dx/2;
%                     total_panels(kk,ll,mm,24)=3;
%                     bndry_z_u(dumzu,1:3)=total_panels(kk,ll,mm,21:23); 
%                     bndry_z_u(dumzu,4)=3;
%                 end
%                 
%                 if ((kk-1)>=1 && boolean_tens(kk-1,ll,mm)<1e-12) || kk==1 %left bndry
%                     dumxl=dumxl+1;
%                     bndry_x_l(dumxl,1:3)=total_panels(kk,ll,mm,1:3);
%                     bndry_x_l(dumxl,4)=1;
%                 end
%                 
%                 if ((ll-1)>=1 && boolean_tens(kk,ll-1,mm)<1e-12)|| ll==1%front bndry
%                     dumyf=dumyf+1;
%                     bndry_y_f(dumyf,1:3)=total_panels(kk,ll,mm,9:11);
%                     bndry_y_f(dumyf,4)=2;
%                 end
%                 if ((mm-1)>=1 && boolean_tens(kk,ll,mm-1)<1e-12) || mm==1 %bottom bndry
%                     dumzb=dumzb+1;
%                     bndry_z_b(dumzb,1:3)=total_panels(kk,ll,mm,17:19);
%                     bndry_z_b(dumzb,4)=3;
%                 end
%             end
%         end  
%     end
% end
bndry_x_r=bndry_x_r(1:dumxr,:);
bndry_x_l=bndry_x_l(1:dumxl,:);
bndry_y_f=bndry_y_f(1:dumyf,:);
bndry_y_b=bndry_y_b(1:dumyb,:);
bndry_z_b=bndry_z_b(1:dumzb,:);
bndry_z_u=bndry_z_u(1:dumzu,:);

bndry_panel=[bndry_x_r; bndry_x_l; bndry_y_b; bndry_y_f; bndry_z_u; bndry_z_b];
% bndry_panel(:,12)=1:size(bndry_panel,1);
% if fl_precond==1
%     boolean_tens_cell=cell(5,1);
%     [boolean_tens_cell{1},boolean_tens_cell{2},boolean_tens_cell{3},boolean_tens_cell{4},boolean_tens_cell{5}] = Tucker(boolean_tens,1e-4);
% else
%     boolean_tens_cell=1;
% end
%%%-------------------------------------------------------------------
%                  obtain numbering for geom panels
%%%-------------------------------------------------------------------
% % generate ids of x_directed panels
% tensor_ijk_xdir_pnl=zeros(L+1,M,N);
% dum_cnt=1;
% for mm=1:N% z-variation
%     for ll=1:M % y-variation
%         for kk=1:L+1 % x-variation
%             tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
%             dum_cnt=dum_cnt+1;
%         end
%     end
% end
% % generate ids of  y_directed panels
% tensor_ijk_ydir_pnl=zeros(L,M+1,N);
% for mm=1:N % z-variation
%     for ll=1:M+1 % y-variation
%         for kk=1:L % x-variation
%             tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
%             dum_cnt=dum_cnt+1;
%         end
%     end
% end
% 
% % generate ids of z_directed panels
% tensor_ijk_zdir_pnl=zeros(L,M,N+1);
% for mm=1:N+1 % z-variation
%     for ll=1:M % y-variation
%         for kk=1:L % x-variation
%             tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
%             dum_cnt=dum_cnt+1;
%         end
%     end
% end
[tensor_ijk_xdir_pnl,tensor_ijk_ydir_pnl,tensor_ijk_zdir_pnl]=panel_numbering(L,M,N);

bbox_origin=[0,0,0];
org_x_panel=[bbox_origin(1) bbox_origin(2)+dx/2 bbox_origin(3)+dx/2];
org_y_panel=[bbox_origin(1)+dx/2 bbox_origin(2) bbox_origin(3)+dx/2];
org_z_panel=[bbox_origin(1)+dx/2 bbox_origin(2)+dx/2 bbox_origin(3)];

% for unique panels


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

% disp('Done... Extracting unique panels of a grid')
% disp('-----------------------------------------------------')