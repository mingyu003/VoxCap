clc
close all
clear all

pre_define_the_path_for_folders

% This routine loads the mat file generated at the end of simulation
% and prints/plots the output data

disp('-----------------------------------------------------')
disp('Post-processing...')

% ------------------------------------------------------------------------
%                     Validation with decreasing dx
% -------------------------------------------------------------------------

set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
unknown=[2.0000000000E+01 4.0000000000E+01 5.0000000000E+01 1.0000000000E+02];

FastCap=[3.5558479337E-02 1.3915070548E-02 9.0445364832E-03 6.8666010336E-03];
VoxCap=[3.3898557371E-02 1.2520493223E-02 7.6840511333E-03 5.5708309151E-03];
% figure('visible','off')
FigHandle = figure%('visible','off');

set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
h=semilogy(unknown,FastCap,'b-o'); set(h,'LineWidth',2);
hold on
h=semilogy(unknown,VoxCap,'r-+'); set(h,'LineWidth',2);

legend('FastCap','VoxCap');
axis tight;grid on;xlabel('\Delta{\itv}^{-1} (1/m)');ylabel('\iterr');
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
ylim([1e-3 1e-1])

print('results_numexA_coated_sphere/vali_dec_dx',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Validation with increasing er
% -------------------------------------------------------------------------
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
unknown=[2 2e1 2e2 2e3 2e4 2e5 2e6 2e7];

FastCap=[5.45949336446059E-02 4.33618721188461E+00 6.34880010697360E+01 7.74785134206549E+02 7.78206972932868E+03 7.78547756446068E+04 5.93062673876664E+05 5.93065104508593E+06];
VoxCap=[3.38986446459156E-02 1.27283921451691E-02 1.03091382147611E-02 1.00706016357349E-02 1.00467974173211E-02 1.00443651699599E-02 1.00435973724522E-02 1.00382773887361E-02];
% figure('visible','off')
FigHandle = figure%('visible','off');

set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
h=loglog(unknown,FastCap,'b-o'); set(h,'LineWidth',2);
hold on
h=loglog(unknown,VoxCap,'r-+'); set(h,'LineWidth',2);
legend('FastCap','VoxCap');
axis tight;grid on;
ylabel('\iterr');
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
ylim([1e-3 6e6])
yticks([10^(-2) 10^2 10^6])
print('results_numexA_coated_sphere/vali_inc_er',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Plotting Charge Distribution
% -------------------------------------------------------------------------

load('results_numexA_coated_sphere/data_solution.mat');
load('results_numexA_coated_sphere/data_geo.mat');
disp('-----------------------------------------------------')


% for large structure, maybe you have to plot part by part. 

plot_charges_on_3D_structure(dx,geom_bndry_panels,inds_glob,real(q_charge_vect))


print('results_numexA_coated_sphere/charge_dist',  '-dpng', '-r300')


disp(['Done... Plotting Charge Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

