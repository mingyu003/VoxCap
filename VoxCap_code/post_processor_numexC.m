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

figure;
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

unknown=1:60;
VoxCap=1e-6*abs([[3.48730104041269e-07,-8.81221359629655e-08,-1.03413350467003e-09,-4.36030924749924e-10,-4.06989710878813e-10, ...
    -2.36127210053391e-07,-1.15313040452984e-08,-4.97735615730717e-11,-2.22464219707089e-11,-6.45125750712269e-11,-3.58358683359228e-09, ...
    -2.29780513687435e-10,-3.09295089436486e-11,-1.81994150936063e-11,-4.77439295871904e-11,-8.60405478164122e-10,-5.03007892247347e-11, ...
    -2.34622012766234e-11,-1.38287129119635e-11,-3.87919576972074e-11,-4.27661072779858e-10,-3.18655909055732e-11,-1.76966384984689e-11, ...
    -1.13719207648959e-11,-3.43641183334478e-11,-2.79603046011795e-10,-2.24795696767650e-11,-1.41176469423262e-11,-9.81091804139460e-12, ...
    -3.09707193208118e-11,-2.02413121401226e-10,-1.59258685143209e-11,-1.22902606313751e-11,-8.80279874624951e-12,-2.85852760578509e-11, ...
    -1.54889967535024e-10,-1.45110920009009e-11,-1.15630029905552e-11,-7.75454744390138e-12,-2.70105584690290e-11,-1.25596732184794e-10, ...
    -1.19054266869255e-11,-8.84579493447283e-12,-7.88568588802596e-12,-2.60500970482275e-11,-1.06993766883047e-10,-1.01025926160437e-11, ...
    -8.64460917155803e-12,-6.52842257894073e-12,-2.68375928325571e-11,-9.97878032093546e-11,-1.13008753636568e-11,-8.44088909174883e-12, ...
    -7.91408575540756e-12,-2.98785404336466e-11,-2.40391620319925e-10,-6.97431008726133e-11,-5.54267205996107e-11,-5.12289528691176e-11, ...
    -1.00877491761748e-10]]); 
Order2=1e-21*abs([3.47823e+08 -8.77578e+07      -759340      -421165      -299198 -2.35646e+08 -1.10542e+07      -460736     -75246.7 ...
    -21416 -4.17248e+06       305090      -261558      -100670      91684.7      -552193       114585      -225181     -39414.2      26803.3 ...
    -1.05251e+06      91762.6     -94818.9     -92238.7     -42602.7       530111      -142935      1047.88      7796.65     -19665.9      ...
    -196730     -39697.4      19131.9     -68664.2     -40897.1      -156693     -29391.7     -5006.55      9357.72     -94085.5     -84270.4 ...    
    -3876.57     -4112.84      6245.49       -56371      -131537     -9079.35      2567.78      9627.18     -58049.4      25814.3   ...
    -119377      -143143     -58768.1     -44121.1      -224209     -74428.7     -62710.3     -63372.2     -99316.4]);
Order4=1e-21*abs([3.46648e+08 -8.76567e+07 -1.06202e+06      -441540      -418865 -2.34935e+08 -1.14602e+07     -72770.2     -33210.7 ...
    -70487.5 -3.57794e+06      -248186      13692.2     -21147.4       -46617      -877023     -56508.3       -26722     -14750.7     -35123.6 ...
    -420808     -29506.4     -14873.3     -8476.74     -32070.7      -287728     -26208.1     -14767.2     -9701.32     -29652.2      -204293 ...
    -18497.9     -13620.7     -9121.88     -28757.3      -159678     -15268.7       -11377     -7976.93     -29121.3      -128267     -13046.8 ...
    -9876.18      -7635.6     -28004.6      -110473     -12848.2     -9048.95     -7285.76     -28069.8     -96501.4     -96501.4     -11472.7 ...
    -10500.6     -30622.6      -251135     -71545.9     -57368.3     -53082.6      -104405]);
Order6=1e-21*abs([3.44967e+08   -8.778e+07 -1.06152e+06      -447403      -419528 -2.34098e+08 -1.13764e+07     -42317.6     -22289.2 ...
    -66973.9 -3.62041e+06      -231195     -30655.4     -16649.2     -48716.7      -870551     -55193.1     -23287.7     -13951.5     -39618.1 ...
    -435187     -32343.6     -18211.9     -11857.5       -34371      -283999     -22278.2     -14478.2     -10393.6     -31391.9      -205253  ...
    -17451.5     -12132.8     -9187.19     -29539.4      -158750     -14270.4     -10454.6     -8288.03     -28206.8      -129049     -12226.7 ...
    -9236.78     -7677.04     -27450.5      -110122       -11071     -8739.69     -7351.24     -27578.2      -105083     -10293.7     -7863.12  ...
    -7959.66     -30443.8      -250980     -71930.3     -57560.2     -53204.4      -104774]);
h=semilogy(unknown,VoxCap,'b-o');
set(h,'LineWidth',2); hold on;

h=semilogy(unknown,Order2,'r--*');
set(h,'LineWidth',2); 

legend('VoxCap','FastCap-2^{nd} order multi. exp')

xlabel('Row ID')
ylabel('Capacitance (F)')
ytickformat('%.0e')
yticks([1e-5 1e-3 1e-1]*1e-12);
grid on
% axis tight
xlim([1 60])
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
set(gca,'yminorgrid','off');
print('results_numexC_parallel_buses/second_order',  '-dpng', '-r300')

figure;
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

h=semilogy(unknown,VoxCap,'b-o');
set(h,'LineWidth',2); hold on;

h=semilogy(unknown,Order4,'r--*');
set(h,'LineWidth',2); hold on;

legend('VoxCap','FastCap-4^{th} order multi. exp')

xlabel('Row ID')
ylabel('Capacitance (F)')
ytickformat('%.0e')
yticks([1e-5 1e-3 1e-1]*1e-12);
grid on
% axis tight
xlim([1 60])
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
set(gca,'yminorgrid','off');
print('results_numexC_parallel_buses/forth_order',  '-dpng', '-r300')

figure;
FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(FigHandle, 'Position', [50, 0, 1280, 1024]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

h=semilogy(unknown,VoxCap,'b-o');
set(h,'LineWidth',2); hold on;

h=semilogy(unknown,Order6,'r--*');
set(h,'LineWidth',2); hold on;

legend('VoxCap','FastCap-6^{th} order multi. exp')

xlabel('Row ID')
ylabel('Capacitance (F)')
ytickformat('%.0e')
yticks([1e-5 1e-3 1e-1]*1e-12);
grid on
% axis tight
xlim([1 60])
set(gca,'FontSize',24)
set(gca,'FontName','Times New Roman')
set(gca,'yminorgrid','off');
print('results_numexC_parallel_buses/sixth_order',  '-dpng', '-r300')
% ------------------------------------------------------------------------
%                     Plotting Charge Distribution
% -------------------------------------------------------------------------

load('results_numexC_parallel_buses/data_solution.mat');
load('results_numexC_parallel_buses/data_geo.mat');
disp('-----------------------------------------------------')


% for large structure, maybe you have to plot part by part. 

plot_charges_on_3D_structure(dx,geom_bndry_panels,inds_glob,real(q_charge_vect))


print('results_numexA_coated_sphere/charge_dist',  '-dpng', '-r300')


disp(['Done... Plotting Charge Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

