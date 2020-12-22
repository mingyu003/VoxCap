function plot_charges_on_3D_structure_large(dx,geom_bndry_panels,inds_glob,q_charge_vect)

% x_aligned_cond_inds_temp=cell(num_cond,1);y_aligned_cond_inds_temp=cell(num_cond,1);z_aligned_cond_inds_temp=cell(num_cond,1);
% for ii=1:num_cond
%     x_aligned_cond_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii);
%     y_aligned_cond_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii);
%     z_aligned_cond_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0 & geom_bndry_panels(:,6) == ii);
% end
% x_aligned_diel_inds_temp=cell(num_diel,1); y_aligned_diel_inds_temp=cell(num_diel,1);  z_aligned_diel_inds_temp=cell(num_diel,1);    
% for ii=1:num_diel
%     x_aligned_diel_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == ii);
%     y_aligned_diel_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == ii);
%     z_aligned_diel_inds_temp{ii}=find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == ii);
% end
% 
% x_aligned_cond_inds=cell2mat(x_aligned_cond_inds_temp);x_aligned_diel_inds=cell2mat(x_aligned_diel_inds_temp);
% y_aligned_cond_inds=cell2mat(y_aligned_cond_inds_temp);y_aligned_diel_inds=cell2mat(y_aligned_diel_inds_temp);
% z_aligned_cond_inds=cell2mat(z_aligned_cond_inds_temp);z_aligned_diel_inds=cell2mat(z_aligned_diel_inds_temp);
% 
% clear x_aligned_cond_inds_temp y_aligned_cond_inds_temp z_aligned_cond_inds_temp x_aligned_diel_inds_temp ...
%     y_aligned_diel_inds_temp z_aligned_diel_inds_temp
% 
% 
% inds_sorted=[x_aligned_cond_inds;x_aligned_diel_inds;y_aligned_cond_inds;y_aligned_diel_inds;z_aligned_cond_inds;z_aligned_diel_inds];
% 
% clear x_aligned_cond_inds x_aligned_diel_inds y_aligned_cond_inds y_aligned_diel_inds z_aligned_cond_inds z_aligned_diel_inds



% options for plotting
fl_no_edge_line = 1; % visualize without voxel edge lines
fl_title_on = 1; % add title to the figure
% The constants related to cut
cut_on = 1; % take a cut on the geometry
cut_axis=1; %1->x,2->y,3->z, cut parallel to axis x,y,or z
cut_pnt=715;%2105; % pnt through which cut passes
cut_below=1; % visualize below the cut [1] (<cut_pnt) or above the cut [0] (>cut_pnt)
% The constants related to log scale
log_on = 1 ;% computed dB after normalizing with maximum current
dB_down = 100;




fl_diel=0;
if (sum(geom_bndry_panels(:,7)) > 1e-13); % check whether we have a dielectric panel
    fl_diel = 1; % yes, we have
end

% 1) Sort the result due to numbering in geom_bndry_panels
% if (fl_diel == 1)
%     inds_sorted=[[inds_glob(1,1):inds_glob(1,2)] [inds_glob(3,1):inds_glob(3,2)] [inds_glob(5,1):inds_glob(5,2)]...
%         [inds_glob(2,1):inds_glob(2,2)] [inds_glob(4,1):inds_glob(4,2)] [inds_glob(6,1):inds_glob(6,2)]];
% else
%     inds_sorted=[[inds_glob(1,1):inds_glob(1,2)] [inds_glob(2,1):inds_glob(2,2)] [inds_glob(3,1):inds_glob(3,2)]];
% end

num_cond=15;
num_diel=11;
num_interface=5;
num_tot=num_cond+num_diel+num_interface;
inds_sorted=[];
for ii=1:num_cond
      inds_sorted=[inds_sorted [inds_glob(ii,1):inds_glob(ii,2)]  [inds_glob(num_tot+ii,1):inds_glob(num_tot+ii,2)]  [inds_glob(2*num_tot+ii,1):inds_glob(2*num_tot+ii,2)]];
end
for ii=1:num_diel
    if inds_glob(num_cond+ii,2)-inds_glob(num_cond+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(num_cond+ii,1):inds_glob(num_cond+ii,2)]];
    end
    if inds_glob(num_tot+num_cond+ii,2)-inds_glob(num_tot+num_cond+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(num_tot+num_cond+ii,1):inds_glob(num_tot+num_cond+ii,2)]];
    end
    if inds_glob(2*num_tot+num_cond+ii,2)-inds_glob(2*num_tot+num_cond+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(2*num_tot+num_cond+ii,1):inds_glob(2*num_tot+num_cond+ii,2)]];
    end
end
for ii=1:num_interface
    if inds_glob(num_cond+num_diel+ii,2)-inds_glob(num_cond+num_diel+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(num_cond+num_diel+ii,1):inds_glob(num_cond+num_diel+ii,2)]];
    end
    if inds_glob(num_tot+num_cond+num_diel+ii,2)-inds_glob(num_tot+num_cond+num_diel+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(num_tot+num_cond+num_diel+ii,1):inds_glob(num_tot+num_cond+num_diel+ii,2)]];
    end
    if inds_glob(2*num_tot+num_cond+num_diel+ii,2)-inds_glob(2*num_tot+num_cond+num_diel+ii,1)>1-1e-13
      inds_sorted=[inds_sorted [inds_glob(2*num_tot+num_cond+num_diel+ii,1):inds_glob(2*num_tot+num_cond+num_diel+ii,2)]];
    end
end
inds_sorted=unique(inds_sorted);

% Note -> charge 
%geom_bndry_panels=geom_bndry_panels(inds_sorted,:);
q_charge_sorted = q_charge_vect(inds_sorted);


geom_bndry_panels(:,8:10)=[];
% aa=find(geom_bndry_panels(:,2)>66 & geom_bndry_panels(:,2)<=450 & geom_bndry_panels(:,3)<=750);
% aa=find(geom_bndry_panels(:,2)>450 & geom_bndry_panels(:,2)<=850 & geom_bndry_panels(:,3)<=750);
% aa=find(geom_bndry_panels(:,2)>850 & geom_bndry_panels(:,2)<=1250 & geom_bndry_panels(:,3)<=750);
% aa=find(geom_bndry_panels(:,2)>1250 & geom_bndry_panels(:,2)<=1650 & geom_bndry_panels(:,3)<=750);
% aa=find(geom_bndry_panels(:,2)>1650 & geom_bndry_panels(:,3)<=750);


aa=find(geom_bndry_panels(:,2)>65 & geom_bndry_panels(:,3)<=750);
geom_bndry_panels=geom_bndry_panels(aa,:);
q_charge_sorted=q_charge_sorted(aa,:);
num_unk=size(geom_bndry_panels,1);


%maxx_charge=max(q_charge_sorted);
%q_charge_norm=q_charge_sorted./maxx_charge;

%  2) create a tensor for vertices and store the contribution from each panel
%  to the the vertices in this panel
clear q_charge_vect inds_glob inds_sorted

bb_org=min(geom_bndry_panels(:,1:3));
%bb_org(2)=bb_org(2)-dx/2;
num_elem_xyz=round(((max(geom_bndry_panels(:,1:3))-bb_org)/dx)+1);

vert_tensor=zeros(num_elem_xyz(1),num_elem_xyz(2),num_elem_xyz(3));
vert_tensor_num_edge=zeros(num_elem_xyz(1),num_elem_xyz(2),num_elem_xyz(3));
% each vertex is shared by how many edges

for kk=1:num_unk
    %[kk num_unk]
    tmp_pnt=geom_bndry_panels(kk,1:4);
    verts=[];
    if (abs(tmp_pnt(4)) == 1) % x-directed
        verts=[tmp_pnt(1) tmp_pnt(2)-dx/2 tmp_pnt(3)-dx/2; ...
            tmp_pnt(1) tmp_pnt(2)+dx/2 tmp_pnt(3)-dx/2; ...
            tmp_pnt(1) tmp_pnt(2)+dx/2 tmp_pnt(3)+dx/2; ...
            tmp_pnt(1) tmp_pnt(2)-dx/2 tmp_pnt(3)+dx/2;];
    elseif (abs(tmp_pnt(4)) == 2) % y-directed
        verts=[tmp_pnt(1)-dx/2 tmp_pnt(2) tmp_pnt(3)-dx/2; ...
            tmp_pnt(1)+dx/2 tmp_pnt(2) tmp_pnt(3)-dx/2; ...
            tmp_pnt(1)+dx/2 tmp_pnt(2) tmp_pnt(3)+dx/2; ...
            tmp_pnt(1)-dx/2 tmp_pnt(2) tmp_pnt(3)+dx/2;];
    elseif (abs(tmp_pnt(4)) == 3) % z-directed
        verts=[tmp_pnt(1)-dx/2 tmp_pnt(2)-dx/2 tmp_pnt(3); ...
            tmp_pnt(1)+dx/2 tmp_pnt(2)-dx/2 tmp_pnt(3); ...
            tmp_pnt(1)+dx/2 tmp_pnt(2)+dx/2 tmp_pnt(3); ...
            tmp_pnt(1)-dx/2 tmp_pnt(2)+dx/2 tmp_pnt(3);];
    end
    
    for ll=1:4
        
        inds_ijk=round((verts(ll,1:3)-bb_org)/dx)+1; % vertex ll
        
        vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3)) = ...
            vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3)) + ...
            q_charge_sorted(kk);
        
        vert_tensor_num_edge(inds_ijk(1),inds_ijk(2),inds_ijk(3)) = ...
            vert_tensor_num_edge(inds_ijk(1),inds_ijk(2),inds_ijk(3)) + 1 ;
    end
    
end
clear verts q_charge_sorted
% finding the average potential on each vertex
for kk=1:num_elem_xyz(1)
    for ll=1:num_elem_xyz(2)
        for mm=1:num_elem_xyz(3)
            if (vert_tensor_num_edge(kk,ll,mm) > 0)
                vert_tensor(kk,ll,mm) = vert_tensor(kk,ll,mm) / ...
                    vert_tensor_num_edge(kk,ll,mm);
            end
        end
    end
end
clear vert_tensor_num_edge
if (log_on == 1) % for normalization, finding the maximum current
    maxx_charge=max(max(max(vert_tensor)));
end


% plotting
figure('visible','off')
FigHandle = figure('visible','off');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
for kk=1:num_unk %1601:1610 %1601:1650 % num_unk
    %[kk num_unk]
    tmp_pnt=geom_bndry_panels(kk,1:4);
    verts=[];
    if (abs(tmp_pnt(4)) == 1) % x-directed
        verts=[tmp_pnt(1) tmp_pnt(2)-dx/2 tmp_pnt(3)-dx/2; ...
            tmp_pnt(1) tmp_pnt(2)+dx/2 tmp_pnt(3)-dx/2; ...
            tmp_pnt(1) tmp_pnt(2)+dx/2 tmp_pnt(3)+dx/2; ...
            tmp_pnt(1) tmp_pnt(2)-dx/2 tmp_pnt(3)+dx/2;];
    elseif (abs(tmp_pnt(4)) == 2) % y-directed
        verts=[tmp_pnt(1)-dx/2 tmp_pnt(2) tmp_pnt(3)-dx/2; ...
            tmp_pnt(1)+dx/2 tmp_pnt(2) tmp_pnt(3)-dx/2; ...
            tmp_pnt(1)+dx/2 tmp_pnt(2) tmp_pnt(3)+dx/2; ...
            tmp_pnt(1)-dx/2 tmp_pnt(2) tmp_pnt(3)+dx/2;];
    elseif (abs(tmp_pnt(4)) == 3) % z-directed
        verts=[tmp_pnt(1)-dx/2 tmp_pnt(2)-dx/2 tmp_pnt(3); ...
            tmp_pnt(1)+dx/2 tmp_pnt(2)-dx/2 tmp_pnt(3); ...
            tmp_pnt(1)+dx/2 tmp_pnt(2)+dx/2 tmp_pnt(3); ...
            tmp_pnt(1)-dx/2 tmp_pnt(2)+dx/2 tmp_pnt(3);];
    end
    
    
    inds_ijk=round((verts(1,1:3)-bb_org)/dx)+1; % vertex 1
    val_v1=vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3));
    
    inds_ijk=round((verts(2,1:3)-bb_org)/dx)+1; % vertex 2
    val_v2=vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3));
    
    inds_ijk=round((verts(3,1:3)-bb_org)/dx)+1; % vertex 3
    val_v3=vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3));
    
    inds_ijk=round((verts(4,1:3)-bb_org)/dx)+1; % vertex 4
    val_v4=vert_tensor(inds_ijk(1),inds_ijk(2),inds_ijk(3));
    
    verts(5,1:3)=tmp_pnt(1:3);
    val_v5=(val_v1+val_v2+val_v3+val_v4)/4;
    all_vals=([val_v1; val_v2; val_v3; val_v4; val_v5;]);
    
    % Attention: The patch command doesn't work and doesn't visualize the
    % vertical squares when the number of patches is high. To fix this problem, each square
    % is divided into 4 triangles and this triangles are fed to 'patch'
    % command. The abovementioned problem doesn't occur when triangles are used
    % instead of squares.Of course, this increases the visualization time.
    
    if (cut_on == 0) % No cut
        
        if (log_on == 1)
            %tmp_vect=[val_v1;val_v2;val_v3;val_v4];
            %tmp_vect=max(20*log10(tmp_vect/maxx_charge),-dB_down);
            %ppp=patch(verts(:,1),verts(:,2),verts(:,3),tmp_vect);
            tmp_vect=max(20*log10(abs(all_vals)/maxx_charge),-dB_down);
            ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),tmp_vect([1 2 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),tmp_vect([2 3 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),tmp_vect([3 4 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),tmp_vect([4 1 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            
        else
            %ppp=patch(verts(:,1),verts(:,2),verts(:,3),[val_v1;val_v2;val_v3;val_v4]);
            %hold on
            ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),all_vals([1 2 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),all_vals([2 3 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),all_vals([3 4 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),all_vals([4 1 5])); hold on
            if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
            
        end
        
        if (fl_no_edge_line == 1)
            set(ppp,'EdgeColor','none')
        end
        
    else % if there is a cut
        if (cut_below==1) % visualize below the cut
            
            if (verts(1,cut_axis) <= cut_pnt && verts(2,cut_axis) <= cut_pnt && ...
                    verts(3,cut_axis) <= cut_pnt && verts(4,cut_axis) <= cut_pnt)
                
                
                if (log_on == 1)
                    %tmp_vect=[val_v1;val_v2;val_v3;val_v4];
                    %tmp_vect=max(20*log10(tmp_vect/maxx_charge),-dB_down);
                    %ppp=patch(verts(:,1),verts(:,2),verts(:,3),tmp_vect);
                    tmp_vect=max(20*log10(abs(all_vals)/maxx_charge),-dB_down);
                    ppp=patch(verts([1 2 3 4],1),verts([1 2 3 4],2),verts([1 2 3 4],3),tmp_vect([1 2 3 4])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
%                     ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),tmp_vect([1 2 5])); hold on
%                     if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
%                     ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),tmp_vect([2 3 5])); hold on
%                     if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
%                     ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),tmp_vect([3 4 5])); hold on
%                     if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
%                     ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),tmp_vect([4 1 5])); hold on
%                     if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                else
                    %ppp=patch(verts(:,1),verts(:,2),verts(:,3),[val_v1;val_v2;val_v3;val_v4]);
                    ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),all_vals([1 2 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),all_vals([2 3 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),all_vals([3 4 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),all_vals([4 1 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                end
                
            end
            
        else % visualize above the cut
            
            if (verts(1,cut_axis) >= cut_pnt && verts(2,cut_axis) >= cut_pnt && ...
                    verts(3,cut_axis) >= cut_pnt && verts(4,cut_axis) >= cut_pnt)
                
                if (log_on == 1)
                    %tmp_vect=[val_v1;val_v2;val_v3;val_v4];
                    %tmp_vect=max(20*log10(tmp_vect/maxx_charge),-dB_down);
                    %ppp=patch(verts(:,1),verts(:,2),verts(:,3),tmp_vect);
                    tmp_vect=max(20*log10(abs(all_vals)/maxx_charge),-dB_down);
                    ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),tmp_vect([1 2 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),tmp_vect([2 3 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),tmp_vect([3 4 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),tmp_vect([4 1 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                else
                    %ppp=patch(verts(:,1),verts(:,2),verts(:,3),[val_v1;val_v2;val_v3;val_v4]);
                    %ppp=trisurf([1 2 3; 3 4 1],verts(:,1),verts(:,2),verts(:,3),[val_v1;val_v2]);
                    %ppp=patch(verts(:,1),verts(:,2),verts(:,3),[val_v1;val_v2;val_v3;val_v4]);
                    ppp=patch(verts([1 2 5],1),verts([1 2 5],2),verts([1 2 5],3),all_vals([1 2 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([2 3 5],1),verts([2 3 5],2),verts([2 3 5],3),all_vals([2 3 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([3 4 5],1),verts([3 4 5],2),verts([3 4 5],3),all_vals([3 4 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                    ppp=patch(verts([4 1 5],1),verts([4 1 5],2),verts([4 1 5],3),all_vals([4 1 5])); hold on
                    if (fl_no_edge_line == 1); set(ppp,'EdgeColor','none'); end
                end
                
                if (fl_no_edge_line == 1)
                    set(ppp,'EdgeColor','none')
                end
                
            end
        end
        
    end
    
end
clear geom_bndry_panels verts vert_tensor
map='hot';
colormap(map);
axis tight; grid on; % this should be here, otherwise it doesn't add vertical patches
axis equal; grid on;
view(45,30)
xlabel('\itx');ylabel('\ity');zlabel('\itz');
%set(gcf,'Renderer','opengl');
%set(gcf,'Renderer','zbuffer')
set(gcf,'Renderer','painters')
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
colorbar
if (log_on == 1)
    colorbar; caxis([-dB_down 0]);
end

%set(ppp,'facealpha',1.0);
%set(ppp,'FaceColor','interp');
if (fl_no_edge_line == 1)
    set(ppp,'EdgeColor','none')
end
if(fl_title_on==1);
    if (log_on == 1)
        title('Normalized Charge (dB)');
    else
        title('Charge (C)');
    end;
end
%daspect([217 217 78])
% daspect([36 80 80])
print(gcf,'-djpeg','-r300','charge.jpg');
