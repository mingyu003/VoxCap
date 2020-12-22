function [int_res] = analytical_derivation_1overR(dx,src_cen,obs_cen,unit_normal_src,unit_normal_obs,epsa,epsb)
%-----------------------------------------------
%analytic formula to calculate 1/R integral
%----------------------------------------------
%%INPUTS
%dx: resolution namely the dimension of every voxel
%num_smpl: the number of integral sample points
%src_cen: center coordinates of source panel
%obs_cen: center coordinates of observer panel
%unit_normal_src: unit normal vector of source panel
%,unit_normal_obs: unit normal vector of observer panel
%%OUTPUT
%int_res: the result of i/R integral
%------------------------------------------------

%clc;clear all;close all;
eps2=1e-37;
eps0=8.854187817e-12;
%eps_dom=1e-37;

fl_self_term = 0;
if (abs(obs_cen(1)-src_cen(1))<eps && abs(obs_cen(2)-src_cen(2))<eps && ...
        abs(obs_cen(3)-src_cen(3))<eps && (abs(unit_normal_src(1))-abs(unit_normal_obs(1)))<eps && ...
        (abs(unit_normal_src(2))-abs(unit_normal_obs(2)))<eps && (abs(unit_normal_src(3))-abs(unit_normal_obs(3)))<eps)
    fl_self_term = 1;
end
if (fl_self_term == 1)
    int_res = (dx^2)*2*pi*((epsa+epsb)/(epsa-epsb));
%     int_res = -(dx^2)*2*pi; % interior neumann problem, normal direction
%     pointeing inward
%     int_res = -(dx^2)*2*pi; % exterior neumann problem, normal direction
%     pointeing outward
%       int_res=0;

end
if (fl_self_term == 0)
    
    % 1) Determine on which plane both source and observer are
    % 1a) Check whether surface normals are normalized
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
    
    % 1c) Finding the plane of observer panel
    
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
        fl_parallel=1;
    else
        fl_parallel=0;
    end
    
    if (fl_parallel == 1)
        a_arr=zeros(1,4);
        b_arr=zeros(1,4);
        if (plane_src==1)
            
            z=abs(obs_cen(3)-src_cen(3));
            aij=abs(obs_cen(1)-src_cen(1));
            bij=abs(obs_cen(2)-src_cen(2));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            b_arr(1)=bij-dx;
            b_arr(2)=bij;
            b_arr(3)=bij+dx;
            b_arr(4)=bij;
            CC_w_sgn=obs_cen(3)-src_cen(3);
        elseif (plane_src==2)
            
            z=abs(obs_cen(2)-src_cen(2));
            aij=abs(obs_cen(3)-src_cen(3));
            bij=abs(obs_cen(1)-src_cen(1));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            b_arr(1)=bij-dx;
            b_arr(2)=bij;
            b_arr(3)=bij+dx;
            b_arr(4)=bij;
            
            CC_w_sgn=obs_cen(2)-src_cen(2);
        elseif (plane_src==3)
            
            z=abs(obs_cen(1)-src_cen(1));
            aij=abs(obs_cen(2)-src_cen(2));
            bij=abs(obs_cen(3)-src_cen(3));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            b_arr(1)=bij-dx;
            b_arr(2)=bij;
            b_arr(3)=bij+dx;
            b_arr(4)=bij;
            
            CC_w_sgn=obs_cen(1)-src_cen(1);
        end
        int_res=0;
        for kk=1:4
            for mm=1:4
                rkm=sqrt(a_arr(kk)^2+b_arr(mm)^2+z^2);
                mult_fact=(-1)^(mm+kk);
                % term 1
                if(abs(a_arr(kk)+rkm)<eps2 || abs(rkm)<eps2 )
                    term1_1=0;
                else
                    term1_1=0.5*(b_arr(mm)^2-z^2)*a_arr(kk)*z/(rkm*(a_arr(kk)+rkm));
                end
                if(abs(a_arr(kk)+rkm)<eps2) % this statement CHECKED-CONFIRMED!
                    term1_2=0;
                else
                    term1_2=a_arr(kk)*z*log(a_arr(kk)+rkm);
                end
                term1 = term1_1 - term1_2;
                
                % term2
                if(abs(b_arr(mm)+rkm)<eps2 || abs(rkm)<eps2 )
                    term2_1 =  0;
                else
                    term2_1 =  0.5*(a_arr(kk)^2-z^2)*b_arr(mm)*z/(rkm*(b_arr(mm)+rkm));
                end
                
                if(abs(b_arr(mm)+rkm)<eps2) % this statement CHECKED-CONFIRMED!
                    term2_2=0;
                else
                    term2_2=b_arr(mm)*z*log(b_arr(mm)+rkm);
                end
                term2 = term2_1 - term2_2;
                % term 3
                if(abs(rkm)<eps2) % this statement CHECKED-CONFIRMED!
                    term3_1 = 0;
                else
                    term3_1 = (b_arr(mm)^2-2*z^2+a_arr(kk)^2)*z/(6*rkm);
                end
                
                term3_2 = (2/3)*z*rkm;
                term3 = term3_1 - term3_2;
                % term 4
                
                if(abs(rkm*z)<eps2) % this statement CHECKED-CONFIRMED!
                    term4_1 = 0;
                else
                    term4_1 = a_arr(kk)*b_arr(mm)*atan(a_arr(kk)*b_arr(mm)/(z*rkm));
                end
                
                %term4_1 = a_arr(kk)*b_arr(mm)*atan(a_arr(kk)*b_arr(mm)/(CC*rhoa))
                
                if(abs(rkm)<eps2) || ( abs((z^2*rkm^2)+(a_arr(kk)^2*b_arr(mm)^2)) < 1e-37) % check this if statement
                    term4_2 = 0;
                else
                    term4_2 = (a_arr(kk)^2*b_arr(mm)^2*z/rkm) * (z^2+rkm^2) / (z^2 * rkm^2 + a_arr(kk)^2 * b_arr(mm)^2);
                end
                
                term4 = term4_1 - term4_2;
                
                int_res=int_res+(mult_fact*(term1+term2-term3-term4));
            end
        end
        if (abs(sum(unit_normal_obs)-sign(CC_w_sgn)) > eps2)
            int_res=-int_res;
%         else
%             int_res=-int_res;
        end  
         
           
    elseif(fl_parallel == 0) % use the formula for othogonal panels
        a_arr=zeros(1,4);
        b_arr=zeros(1,2);
        c_arr=zeros(1,2);
        if  plane_src == 1 && plane_obs == 2   % xy,xz
            aij=abs(obs_cen(1)-src_cen(1));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(3)-src_cen(3));
            cij=abs(obs_cen(2)-src_cen(2));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(2)-src_cen(2);
        elseif plane_src == 2 && plane_obs == 1 %xz,xy
            aij=abs(obs_cen(1)-src_cen(1));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(2)-src_cen(2));
            cij=abs(obs_cen(3)-src_cen(3));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(3)-src_cen(3);
        elseif (plane_src == 1 && plane_obs == 3 ) % xy,yz
            aij=abs(obs_cen(2)-src_cen(2));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(3)-src_cen(3));
            cij=abs(obs_cen(1)-src_cen(1));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(1)-src_cen(1);
        elseif ( plane_src == 3 && plane_obs == 1) % yz ,xy
            aij=abs(obs_cen(2)-src_cen(2));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(1)-src_cen(1));
            cij=abs(obs_cen(3)-src_cen(3));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(3)-src_cen(3);
        elseif (plane_src == 2 && plane_obs == 3)% xz,yz
            aij=abs(obs_cen(3)-src_cen(3));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(2)-src_cen(2));
            cij=abs(obs_cen(1)-src_cen(1));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(1)-src_cen(1);
        elseif (plane_src == 3 && plane_obs == 2)% yz, xz
            aij=abs(obs_cen(3)-src_cen(3));
            a_arr(1)=aij-dx;
            a_arr(2)=aij;
            a_arr(3)=aij+dx;
            a_arr(4)=aij;
            bij=abs(obs_cen(1)-src_cen(1));
            cij=abs(obs_cen(2)-src_cen(2));
            b_arr(1)=bij+dx/2;
            b_arr(2)=bij-dx/2;
            c_arr(1)=cij+dx/2;
            c_arr(2)=cij-dx/2;
            CC_w_sgn=obs_cen(2)-src_cen(2);
        end
        int_res=0;
        for kk=1:4
            for mm=1:2
                for ll=1:2
                    rkml=sqrt(a_arr(kk)^2+b_arr(mm)^2+c_arr(ll)^2);
                    mult_fact=(-1)^(mm+kk+ll+1);
                    
                    % term1
                    
                    if(abs(b_arr(mm)+rkml)<eps2 || abs(rkml) < eps2)
                        term1 = 0;
                    else
                        term1_a = 3*(a_arr(kk)^2)*(c_arr(ll)^2)-(c_arr(ll)^4);
                        term1_b = 3*(a_arr(kk)^2-c_arr(ll)^2)*(a_arr(kk)^2+b_arr(mm)^2+c_arr(ll)^2+b_arr(mm)*rkml)*log(b_arr(mm)+rkml);
                        term1_c = 6*rkml*(b_arr(mm)+rkml);
                        term1 = (term1_a + term1_b)/term1_c;
                    end
                    
                    % term 2
                    
                    if(abs(rkml)<eps2)
                        term2 = 0;
                    else
                        term2 = b_arr(mm)*(3*(a_arr(kk)^2)-(b_arr(mm)^2))/(6*rkml);
                    end
                    
                    % term 3
                    
                    if(abs(a_arr(kk)+rkml)<eps2 || abs(rkml)<eps2 )
                        term3_a = 0;
                    else
                        term3_a = (a_arr(kk)*b_arr(mm)*(c_arr(ll)^2))/(rkml*(a_arr(kk)+rkml));
                    end
                    
                    if(abs(a_arr(kk)+rkml)<eps2)
                        term3_b = 0 ;
                    else
                        term3_b = a_arr(kk)*b_arr(mm)*log(a_arr(kk)+rkml);
                    end
                    
                    term3 = term3_a + term3_b;
                    
                    % term 4
                    
                    if(abs(rkml)<eps2)
                        term4 = 0 ;
                    else
                        term4 = b_arr(mm)*(a_arr(kk)^2+b_arr(mm)^2+2*(c_arr(ll)^2)) / (3*rkml);
                    end
                    
                    if(abs(rkml)<eps2 || abs(a_arr(kk)^2+c_arr(ll)^2)<eps2)
                        term5 = 0 ;
                    else
                        term5 = (a_arr(kk)^4*b_arr(mm)) / (6*(a_arr(kk)^2+c_arr(ll)^2)*rkml) ;
                    end
                    
                    if(abs(rkml)<eps2 || abs(b_arr(mm)^2+c_arr(ll)^2)<eps2)
                        term6 = 0;
                    else
                        term6 = (a_arr(kk)^2*b_arr(mm)^3)/(2*(b_arr(mm)^2+c_arr(ll)^2)*rkml);
                    end
                    
                    term7_a = (a_arr(kk)*c_arr(ll)/2);
                    
                    if(abs(rkml)<eps2 || abs(b_arr(mm)^2+c_arr(ll)^2)<eps2 || abs(a_arr(kk)^2+c_arr(ll)^2)<eps2)
                        term7_b = 0;
                    else
                        term7_b1 = (a_arr(kk)*b_arr(mm)*c_arr(ll)*(a_arr(kk)^2+b_arr(mm)^2+2*c_arr(ll)^2));
                        term7_b2 = (a_arr(kk)^2+c_arr(ll)^2)*(b_arr(mm)^2+c_arr(ll)^2)*rkml;
                        term7_b = term7_b1/term7_b2;
                    end
                    
                    if(abs(rkml*c_arr(ll))<eps2)
                        term7_c = 0;
                    else
                        term7_c = 2*atan(a_arr(kk)*b_arr(mm)/(c_arr(ll)*rkml));
                    end
                    
                    term7 = term7_a*(-term7_b+term7_c);
                    
                    int_res=int_res+(mult_fact*(term1+term2+term3-term4-term5-term6-term7));
                end
            end
        end
        if (abs(sum(unit_normal_obs)-sign(CC_w_sgn)) > eps2)
            int_res=-int_res;
%         else
%             int_res=-int_res;
        end
    end
end






