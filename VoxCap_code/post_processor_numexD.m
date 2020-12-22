clc
close all
clear all

pre_define_the_path_for_folders

% This routine loads the mat file generated at the end of simulation
% and prints/plots the output data

disp('-----------------------------------------------------')
disp('Post-processing...')

% ------------------------------------------------------------------------
%                    Comparison with FastCap
% -------------------------------------------------------------------------

unknown=1:45;

FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

vox=1e-9*[[1.50106344231808e-07;-8.35053629636467e-08;-3.16223904138976e-09;-1.08489832954056e-09;-1.04785695782535e-09;-7.70331205787896e-09;-4.41183505914156e-09;-4.31718290432660e-09;-4.41301809154089e-09;-7.79937727778397e-09;-3.50067251049801e-09;-1.19025932868592e-09;-6.46288297107929e-10;-3.86795440286204e-10;-5.26884770590830e-10]];
vox=abs(vox);
fast2=1e-9*1e-9*abs([149.685      -83.869     -3.20879     -1.11541     -1.04719     -11.4333     -6.69084     -6.55147     -6.68934     -11.4716     -3.78663      -1.3621    -0.724592     -0.41085    -0.524701]);
h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;
h=semilogy(unknown,fast2,'r->');
set(h,'LineWidth',2); hold on;
xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
axis tight
legend('VoxCap','FastCap-2^{nd} order multi. exp')
% legend('FastCap-2^{nd} order multi. exp')
ytickformat('%.0e')
% yticks([1e-10 1e-9 1e-8]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
print('results_numexD_multilayer/second_order',  '-dpng', '-r300')

FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');


fast4=1e-9*1e-12*abs([222527      -130047     -3341.62     -1086.22     -618.223     -412.092     -300.311     -236.186     -188.491  ...
    -155.427     -132.786     -116.876     -107.383     -106.638     -265.412     -5275.62     -4384.71     -4419.65     -4431.19  ...
    -4438.36     -4440.95     -4444.04     -4440.48     -4444.04     -4441.21     -4439.08     -4430.81     -4419.64     -4384.84  ...
    -5276.55     -3305.96     -461.948     -265.533     -188.849       -152.8     -126.863     -109.529     -93.9355     -83.9413   ...
    -78.3219     -73.0325     -70.1066     -70.1006     -77.8558     -210.056]);

h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;
h=semilogy(unknown,fast4,'r->');
set(h,'LineWidth',2); hold on;
xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
axis tight
legend('VoxCap','FastCap-4^{th} order multi. exp')
% legend('FastCap-2^{nd} order multi. exp')
ytickformat('%.0e')
% yticks([1e-10 1e-9 1e-8]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
print('results_numexD_multilayer/forth_order',  '-dpng', '-r300')

% return
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
fast6=1e-21*abs([222580      -130083     -3342.67     -1097.29     -619.391     -413.081     -301.347     -233.301       -187.2   ...
    -155.265     -132.887     -117.228     -107.249      -106.59     -265.405     -5276.94     -4386.45     -4420.15     -4432.66  ...
    -4438.86     -4442.28     -4444.06     -4444.28     -4443.95      -4442.3     -4438.99     -4432.42     -4419.18     -4385.22  ...
    -5275.96     -3323.39      -460.01     -262.432     -192.122     -151.515     -126.178     -108.357     -94.3228     -85.2044   ...
    -78.151     -73.1351     -70.0311      -69.919     -77.8047     -210.338]);

h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;
h=semilogy(unknown,fast6,'r->');
set(h,'LineWidth',2); hold on;
xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
axis tight
legend('VoxCap','FastCap-6^{th} order multi. exp')
% legend('FastCap-2^{nd} order multi. exp')
ytickformat('%.0e')
% yticks([1e-10 1e-9 1e-8]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
print('results_numexD_multilayer/sixth_order',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Plotting Charge Distribution
% -------------------------------------------------------------------------

load('results_numexD_multilayer/data_solution.mat');
load('results_numexD_multilayer/data_geo.mat');
disp('-----------------------------------------------------')


% for large structure, maybe you have to plot part by part. 

plot_charges_on_3D_structure_large(dx,geom_bndry_panels,inds_glob,real(q_charge_vect))


print('results_numexD_multilayer/charge_dist',  '-dpng', '-r300')


disp(['Done... Plotting Charge Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

