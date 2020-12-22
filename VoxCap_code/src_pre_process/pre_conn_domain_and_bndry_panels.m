function [num_panels,ids_xyz_aligned]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels)

% Numbers of x, y, anz z aligned panels on computational domain

num_x_aligned_domain = length(find(abs(comp_dom_panels(:,4)) == 1));
num_y_aligned_domain = length(find(abs(comp_dom_panels(:,4)) == 2));
num_z_aligned_domain = length(find(abs(comp_dom_panels(:,4)) == 3));
num_panels_domain = size(comp_dom_panels,1);

% Numbers of x, y, and z aligned panels on geometry

num_x_aligned = length(find(abs(geom_bndry_panels(:,4)) == 1));
num_y_aligned = length(find(abs(geom_bndry_panels(:,4)) == 2));
num_z_aligned = length(find(abs(geom_bndry_panels(:,4)) == 3));

num_panels = size(geom_bndry_panels,1);

% Global IDs of x,y,z aligned geometry panels (due to comp domain panel numbering)

ids_x_aligned_geom = geom_bndry_panels(1:num_x_aligned,5);
ids_y_aligned_geom = geom_bndry_panels(num_x_aligned+1:num_x_aligned+num_y_aligned,5);
ids_z_aligned_geom = geom_bndry_panels(num_x_aligned+num_y_aligned+1:num_x_aligned+num_y_aligned+num_z_aligned,5);


% Local IDs of x,y,z aligned panels (due to comp domain panel numbering)
% among all x,y,z aligned panels of comp domain

ids_xyz_aligned=cell(6,1);

ids_xyz_aligned{1} = ids_x_aligned_geom;
ids_xyz_aligned{2} = ids_y_aligned_geom - num_x_aligned_domain;
ids_xyz_aligned{3} = ids_z_aligned_geom - (num_x_aligned_domain + num_y_aligned_domain);

ids_xyz_aligned{4} = sign(geom_bndry_panels(1:num_x_aligned,4));
ids_xyz_aligned{5} = sign(geom_bndry_panels(num_x_aligned+1:num_x_aligned+num_y_aligned,4));
ids_xyz_aligned{6} = sign(geom_bndry_panels(num_x_aligned+num_y_aligned+1:num_x_aligned+num_y_aligned+num_z_aligned,4));


disp(['# of panels on boundaries ::: ',num2str(num_panels)])
disp(['# of panels on computational domain ::: ',num2str(num_panels_domain)])