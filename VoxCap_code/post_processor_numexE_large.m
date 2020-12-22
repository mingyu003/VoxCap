clc
close all
clear all

pre_define_the_path_for_folders

% This routine loads the mat file generated at the end of simulation
% and prints/plots the output data

disp('-----------------------------------------------------')
disp('Post-processing...')
load('results_numexD_multilayer/data_solution.mat');
load('results_numexD_multilayer/data_geo.mat');
disp('-----------------------------------------------------')
% ------------------------------------------------------------------------
%                    Comparison with FastCap
% -------------------------------------------------------------------------
figure;
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
load data_solution;
unknown=1:100;
VoxCap=1e-3*abs(C_mat);
h=semilogy(unknown,VoxCap,'b-o');
set(h,'LineWidth',2); hold on;

% legend('VoxCap')
xlabel('Row ID')
ylabel('Capacitance (F)')
ytickformat('%.0e')
yticks([1e-14 1e-13 1e-12 1e-11 1e-10]);
grid on
axis tight
xlim([1 100])
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')

print('results_numexE_large/forth_order',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Plotting Charge Distribution
% -------------------------------------------------------------------------


% for large structure, maybe you have to plot part by part. 

plot_charges_on_3D_structure(dx,geom_bndry_panels,inds_glob,real(q_charge_vect))


print('results_numexE_large/charge_dist',  '-dpng', '-r300')


disp(['Done... Plotting Charge Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

