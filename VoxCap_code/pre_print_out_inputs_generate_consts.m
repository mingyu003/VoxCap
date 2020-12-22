
disp('------------------------------------------------------------------')
disp('VoxCap: Capacitance Extraction Simulator for Voxelized Geometries ')
disp('                                                                  ')
disp('      by Mingyu Wang, Abdulkadir C. Yucel, Cheng Qian,            ')
disp('           Nanyang Technological University (NTU)                 ')
disp('                      Jacob K. White                              ')
disp('         Massachusetts Institute of Technology (MIT)              ')
disp('------------------------------------------------------------------')
disp('Inputs for simulation :::')
disp(['Voxel Size = ', num2str(Res)])
disp(['Tolerance for iterative solver = ', num2str(tol)])
disp(['# of maximum inner / outer GMRES iterations = ', num2str(inner_it),' / ',num2str(outer_it)])

disp(['# of seperate conductors = ', num2str(num_cond)])
disp(['# of seperate dielectric = ', num2str(num_diel)])

if (num_cond == 0)
    disp('num_cond should be more than 0 ...')
    error('The VoxCap can not be executed without defining a conductor')
else
    disp('Conductors reside in the media with following permittivities:')
    for kk=1:num_cond
        disp(['Conductor # ',num2str(kk),' is in medium with permittivity ', num2str(Eps_inout(kk,2))])
    end
end

disp(['Attention: The current code will execute the following ::: '])
disp(['1) Single/multiple conductor(s)'])
disp(['2) Single/multiple conductor(s) with single dielectric'])
disp(['The following will be implemented later on!!!'])
disp(['---Single/multiple conductor(s) with multiple dielectrics---'])
disp(['Attention: The Tucker part refered to htucker toolbox developed by C. Tobler and D. Kressner, EPF Lausanne '])

