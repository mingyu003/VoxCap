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
unknown=1:20;

FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

vox=[1.0753130023724E-11 -1.0349104862948E-11 -6.3152853534593E-14 -3.2585004586359E-14 -2.1036240331652E-14 -1.5216435890659E-14 ... 
    -1.1721691333588E-14 -9.4479458616743E-15 -7.8394276438388E-15 -6.7032643050063E-15 -5.8247810799110E-15 -5.1605103357900E-15 ...
-4.6538384413829E-15 -4.2328770547276E-15 -3.9238753851915E-15 -3.6899254261575E-15 -3.5347580368106E-15 -3.5032932465886E-15 ...
-3.6861838902713E-15 -1.5359227109813E-14];
vox=abs(vox);
Fast_2=[1.0694900000000E-11 -1.0307100000000E-11 -4.4285900000000E-14 -3.7012700000000E-14 -2.3094400000000E-14 -1.0586800000000E-14 ...
-1.2146600000000E-14 -9.2099100000000E-15 -7.4917200000000E-15 -6.1285100000000E-15 -6.6330400000000E-15 -5.0760700000000E-15 ...
-4.6645800000000E-15 -4.2733400000000E-15 -3.9591000000000E-15 -3.7811700000000E-15 -3.4277100000000E-15 -3.5013400000000E-15 ...
-3.7698300000000E-15 -1.5443000000000E-14];

Fast_2=abs(Fast_2);

h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;


h=semilogy(unknown,Fast_2,'r--*');
set(h,'LineWidth',2); hold on;

xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
axis tight
legend('VoxCap','FastCap-2^{nd} order multi. exp')
% xticks([1 3 5 7 9 11 13 15 17 19]);
ytickformat('%.0e')
yticks([1e-15 1e-13 1e-11]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
print('results_numexE_small/second_order',  '-dpng', '-r300')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

vox=[1.0753130023724E-11 -1.0349104862948E-11 -6.3152853534593E-14 -3.2585004586359E-14 -2.1036240331652E-14 -1.5216435890659E-14 ... 
    -1.1721691333588E-14 -9.4479458616743E-15 -7.8394276438388E-15 -6.7032643050063E-15 -5.8247810799110E-15 -5.1605103357900E-15 ...
-4.6538384413829E-15 -4.2328770547276E-15 -3.9238753851915E-15 -3.6899254261575E-15 -3.5347580368106E-15 -3.5032932465886E-15 ...
-3.6861838902713E-15 -1.5359227109813E-14];
vox=abs(vox);
Fast_6=[1.06401E-11 -1.02366E-11 -6.28412E-14 -3.2514E-14 -2.10645E-14 -1.52089E-14 -1.17165E-14 -9.43412E-15 -7.84866E-15 -6.68545E-15 ...
-5.8369E-15 -5.15883E-15 -4.63998E-15 -4.23476E-15 -3.91958E-15 -3.68781E-15 -3.53868E-15 -3.49893E-15 -3.67827E-15 -1.5348E-14];

Fast_6=abs(Fast_6);

h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;


h=semilogy(unknown,Fast_6,'r--*');
set(h,'LineWidth',2); hold on;

xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
% axis tight
legend('VoxCap','FastCap-6^{th} order multi. exp')
ytickformat('%.0e')
yticks([1e-15 1e-13 1e-11]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
print('results_numexE_small/sixth_order',  '-dpng', '-r300')

FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

vox=[1.0753130023724E-11 -1.0349104862948E-11 -6.3152853534593E-14 -3.2585004586359E-14 -2.1036240331652E-14 -1.5216435890659E-14 ... 
    -1.1721691333588E-14 -9.4479458616743E-15 -7.8394276438388E-15 -6.7032643050063E-15 -5.8247810799110E-15 -5.1605103357900E-15 ...
-4.6538384413829E-15 -4.2328770547276E-15 -3.9238753851915E-15 -3.6899254261575E-15 -3.5347580368106E-15 -3.5032932465886E-15 ...
-3.6861838902713E-15 -1.5359227109813E-14];
vox=abs(vox);
Fast_4=[1.06426E-11 -1.02419E-11 -6.01419E-14 -3.2764E-14 -2.05829E-14 -1.52429E-14 -1.17535E-14 -9.46399E-15 -7.84552E-15 -6.85455E-15 ...
-5.67962E-15 -5.19063E-15 -4.65228E-15 -4.23191E-15 -3.9413E-15 -3.6647E-15 -3.54535E-15 -3.50465E-15 -3.67616E-15 -1.53424E-14];

Fast_4=abs(Fast_4);

h=semilogy(unknown,vox,'b-o');
set(h,'LineWidth',2); hold on;

h=semilogy(unknown,Fast_4,'r--*');
set(h,'LineWidth',2); hold on;

xlabel('Row ID')
% xlabel('N')
ylabel({'Capacitance (F)'})
grid on
% axis tight
legend('VoxCap','FastCap-4^{th} order multi. exp')
ytickformat('%.0e')
yticks([1e-15 1e-13 1e-11]);
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')

print('results_numexE_small/forth_order',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Plotting Charge Distribution
% -------------------------------------------------------------------------

load('results_numexD_multilayer/data_solution.mat');
load('results_numexD_multilayer/data_geo.mat');
disp('-----------------------------------------------------')


% for large structure, maybe you have to plot part by part. 

plot_charges_on_3D_structure(dx,geom_bndry_panels,inds_glob,real(q_charge_vect))


print('results_numexE_small/charge_dist',  '-dpng', '-r300')


disp(['Done... Plotting Charge Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

