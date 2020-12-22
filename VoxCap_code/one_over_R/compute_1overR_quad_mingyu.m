function [int_res] = compute_1overR_quad_mingyu(dx,src_cen,obs_cen,unit_normal_src,unit_normal_obs,smpl_pnts,smpl_wghts)
%---------------------------------------------------------
%generate 1/R integral matrix with Gauss-Legendre quadrature
%--------------------------------------------------------
%%INPUTS
%dx: resolution namely the dimension of every voxel
%num_smpl: the number of integral sample points
%src_cen: center coordinates of source panel
%obs_cen: center coordinates of observer panel
%unit_normal_src: unit normal vector of source panel
%unit_normal_obs: unit normal vector of observer panel
%smpl_pnts,smpl_wghts: gauss points and weights
%%OUTPUT
%int_res: the result of 1/R integral
%--------------------------------------------------------

eps2=1e-13;
if ((norm(unit_normal_src)-1) > eps2)
    unit_normal_src=unit_normal_src/(norm(unit_normal_src));
end
if ((norm(unit_normal_obs)-1) > eps2)
    unit_normal_obs=unit_normal_obs/(norm(unit_normal_obs));
end
plane_src=0;
if (abs(abs(unit_normal_src(3))-1) < eps2 && abs(unit_normal_src(2)) < eps2 && abs(unit_normal_src(1)) < eps2)
    % on xy plane
    %disp('Source panel is parallel to xy plane');
    plane_src=1;
elseif (abs(abs(unit_normal_src(2))-1) < eps2 && abs(unit_normal_src(3)) < eps2 && abs(unit_normal_src(1)) < eps2)
    % on xz plane
    %disp('Source panel is parallel to xz plane');
    plane_src=2;
elseif (abs(abs(unit_normal_src(1))-1) < eps2 && abs(unit_normal_src(2)) < eps2 && abs(unit_normal_src(3)) < eps2)
    % on yz plane
    %disp('Source panel is parallel to yz plane');
    plane_src=3;
end

if (plane_src==0)
    error('Error! Source panel should be parallel to xy, xz, or yz planes')
end
plane_obs=0;
if (abs(abs(unit_normal_obs(3))-1) < eps2 && abs(unit_normal_obs(2)) < eps2 && abs(unit_normal_obs(1)) < eps2)
    % on xy plane
    %disp('Observer panel is parallel to xy plane');
    plane_obs=1;
elseif (abs(abs(unit_normal_obs(2))-1) < eps2 && abs(unit_normal_obs(3)) < eps2 && abs(unit_normal_obs(1)) < eps2)
    % on xz plane
    %disp('Observer panel is parallel to xz plane');
    plane_obs=2;
elseif (abs(abs(unit_normal_obs(1))-1) < eps2 && abs(unit_normal_obs(2)) < eps2 && abs(unit_normal_obs(3)) < eps2)
    % on yz plane
    %disp('Observer panel is parallel to yz plane');
    plane_obs=3;
end

if (plane_obs==0)
    error('Error! Observer panel should be parallel to xy, xz, or yz planes')
end
if (plane_obs == plane_src)
    %disp('the panels are on the same plane')
    form_type=1;
else
    %disp('the planes of panels are orthogonal')
    form_type=2;
end

%     intervals=[obs_cen(1)-dx/2 obs_cen(1)+dx/2; src_cen(1)-dx/2 src_cen(1)+dx/2; ...
%         obs_cen(2)-dx/2 obs_cen(2)+dx/2; src_cen(2)-dx/2 src_cen(2)+dx/2; ...
%         obs_cen(3)-dx/2 obs_cen(3)+dx/2; src_cen(3)-dx/2 src_cen(3)+dx/2;];
%     int_shift(:) = (intervals(:,2)+intervals(:,1))*0.5;

if (form_type == 1)
    if (plane_src == 1)
        z1=src_cen(3);
        z2=obs_cen(3);
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(1));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(2));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(2));
        Grfn=@(x,xpr,y,ypr) 1./(sqrt(((x-xpr).^2+(y-ypr).^2)+(z2-z1)^2));
%          Grfn=@(x,xpr,y,ypr) 1./(sqrt(((x-xpr).^2+(y-ypr).^2+(intervals(5,1)-intervals(6,1)).^2)));
%          temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(2),smpl_pnts(:,3)+int_shift(3),smpl_pnts(:,4)+int_shift(4));
%          int_res=sum(temp.*smpl_wghts);
    elseif (plane_src == 2)
        y1=src_cen(2);
        y2=obs_cen(2);
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(1));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(3));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
        Grfn=@(x,xpr,z,zpr) 1./(sqrt(((x-xpr).^2+(z-zpr).^2)+(y2-y1)^2));
%             Grfn=@(x,xpr,z,zpr) 1./(sqrt(((x-xpr).^2+(z-zpr).^2+(intervals(4,1)-intervals(3,1)).^2)));
%             temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(2),smpl_pnts(:,3)+int_shift(5),smpl_pnts(:,4)+int_shift(6));
%             int_res=sum(temp.*smpl_wghts);
    elseif (plane_src == 3)
        x1=src_cen(1);
        x2=obs_cen(1);
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(2));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(2));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(3));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
        Grfn=@(y,ypr,z,zpr) 1./(sqrt(((y-ypr).^2+(z-zpr).^2)+(x2-x1)^2));
%          Grfn=@(y,ypr,z,zpr) 1./(sqrt(((y-ypr).^2+(z-zpr).^2+(intervals(2,1)-intervals(1,1)).^2)));  
%          temp=Grfn(smpl_pnts(:,1)+int_shift(3),smpl_pnts(:,2)+int_shift(4),smpl_pnts(:,3)+int_shift(5),smpl_pnts(:,4)+int_shift(6));
%          int_res=sum(temp.*smpl_wghts);  
        
    end
    
elseif(form_type == 2)
    if (plane_src == 1 && plane_obs == 2)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(1));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(2));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(3));
%         Grfn=@(x,xpr,ypr,z) 1./(sqrt(((x-xpr).^2+(obs_cen(2)-ypr).^2)+(z-src_cen(3)).^2));
        Grfn=@(x,xpr,z,ypr) 1./(sqrt(((x-xpr).^2+(obs_cen(2)-ypr).^2)+(z-src_cen(3)).^2));
%         Grfn=@(x,xpr,z,ypr) 1./(sqrt(((x-xpr).^2+(obs_cen(2)-ypr).^2+(z-src_cen(3)).^2)));
%         temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(2),smpl_pnts(:,3)+int_shift(5),smpl_pnts(:,4)+int_shift(4));
%         int_res=sum(temp.*smpl_wghts);    
    elseif (plane_src == 1 && plane_obs == 3)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+src_cen(1));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(2));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(2));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+obs_cen(3));
%         Grfn=@(xpr,y,ypr,z) 1./(sqrt(((obs_cen(1)-xpr).^2+(y-ypr).^2)+(z-src_cen(3)).^2));
        Grfn=@(xpr,ypr,y,z) 1./(sqrt(((obs_cen(1)-xpr).^2+(y-ypr).^2)+(z-src_cen(3)).^2));
%         Grfn=@(xpr,ypr,y,z) 1./(sqrt(((obs_cen(1)-xpr).^2+(y-ypr).^2+(z-src_cen(3)).^2)));
%         temp=Grfn(smpl_pnts(:,1)+int_shift(2),smpl_pnts(:,2)+int_shift(4),smpl_pnts(:,3)+int_shift(3),smpl_pnts(:,4)+int_shift(5));
%         int_res=sum(temp.*smpl_wghts);
    
    elseif (plane_src == 2 && plane_obs == 3)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+src_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+obs_cen(2));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(3));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
        Grfn=@(xpr,y,z,zpr) 1./(sqrt(((obs_cen(1)-xpr).^2+(y-src_cen(2)).^2)+(z-zpr).^2));
%         Grfn=@(xpr,y,z,zpr) 1./(sqrt(((obs_cen(1)-xpr).^2+(y-src_cen(2)).^2+(z-zpr).^2)));
%         temp=Grfn(smpl_pnts(:,1)+int_shift(2),smpl_pnts(:,2)+int_shift(3),smpl_pnts(:,3)+int_shift(5),smpl_pnts(:,4)+int_shift(6));
%         int_res=sum(temp.*smpl_wghts);
    
    elseif (plane_src == 2 && plane_obs == 1)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+src_cen(1));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+obs_cen(2));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
        Grfn=@(x,xpr,y,zpr) 1./(sqrt(((x-xpr).^2+(y-src_cen(2)).^2)+(obs_cen(3)-zpr).^2));
    
%         Grfn=@(x,xpr,y,zpr) 1./(sqrt(((x-xpr).^2+(y-src_cen(2)).^2+(obs_cen(3)-zpr).^2)));
%         temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(2),smpl_pnts(:,3)+int_shift(3),smpl_pnts(:,4)+int_shift(6));
%         int_res=sum(temp.*smpl_wghts);
   
    elseif (plane_src == 3 && plane_obs == 1)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+obs_cen(2));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+src_cen(2));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
        Grfn=@(x,y,ypr,zpr) 1./(sqrt(((x-src_cen(1)).^2+(y-ypr).^2)+(obs_cen(3)-zpr).^2));
%           Grfn=@(x,y,ypr,zpr) 1./(sqrt(((x-src_cen(1)).^2+(y-ypr).^2+(obs_cen(3)-zpr).^2)));
%           temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(3),smpl_pnts(:,3)+int_shift(4),smpl_pnts(:,4)+int_shift(6));
%           int_res=sum(temp.*smpl_wghts);   
    
    elseif (plane_src == 3 && plane_obs == 2)
        smpl_pnts(:,1)=(smpl_pnts(:,1)+obs_cen(1));
        smpl_pnts(:,3)=(smpl_pnts(:,3)+src_cen(2));
        smpl_pnts(:,2)=(smpl_pnts(:,2)+obs_cen(3));
        smpl_pnts(:,4)=(smpl_pnts(:,4)+src_cen(3));
%         Grfn=@(x,ypr,z,zpr) 1./(sqrt(((x-src_cen(1)).^2+(obs_cen(2)-ypr).^2)+(z-zpr).^2));
        Grfn=@(x,z,ypr,zpr) 1./(sqrt(((x-src_cen(1)).^2+(obs_cen(2)-ypr).^2)+(z-zpr).^2));
%           Grfn=@(x,z,ypr,zpr) 1./(sqrt(((x-src_cen(1)).^2+(obs_cen(2)-ypr).^2+(z-zpr).^2)));
%           temp=Grfn(smpl_pnts(:,1)+int_shift(1),smpl_pnts(:,2)+int_shift(5),smpl_pnts(:,3)+int_shift(4),smpl_pnts(:,4)+int_shift(6));
%           int_res=sum(temp.*smpl_wghts);         

    end
end
temp=Grfn(smpl_pnts(:,1),smpl_pnts(:,2),smpl_pnts(:,3),smpl_pnts(:,4));
% int_res=sum(smpl_wghts'*temp);
int_res=sum(temp.*smpl_wghts);

