function [int_res] = compute_1overR_analy_mingyu(dx,src_cen,obs_cen,unit_normal_src,unit_normal_obs)
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
%eps_dom=1e-37;
fl_self_term = 0;
if (abs(obs_cen(1)-src_cen(1))<eps && abs(obs_cen(2)-src_cen(2))<eps && ...
        abs(obs_cen(3)-src_cen(3))<eps && abs(unit_normal_src(1)-unit_normal_obs(1))<eps && ...
        abs(unit_normal_src(2)-unit_normal_obs(2))<eps && abs(unit_normal_src(3)-unit_normal_obs(3))<eps)
    fl_self_term = 1;
end
%when self term, use simplified formula
if (fl_self_term == 1)
    
    const_fafbsasb=(2/3);
    arg_log=(dx+sqrt(2*dx^2))/dx;
    term1=3*dx^3*log(arg_log);
    term2=(2*dx^2)^(3/2);
    term3=term1;
    term4=(2*dx^3);
    
    int_res=const_fafbsasb*(term1-term2+term3+term4);
    
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
    
    if (plane_obs==plane_src)
        fl_parallel=1;
    else
        fl_parallel=0;
    end

    if (fl_parallel == 1)
        if (plane_src==1) 

            z=abs(obs_cen(3)-src_cen(3));
            aij=abs(obs_cen(1)-src_cen(1));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(2)-src_cen(2));
            b_arr=[bij-dx bij bij+dx bij];
            
        elseif (plane_src==2)
            
            z=abs(obs_cen(2)-src_cen(2));
            aij=abs(obs_cen(3)-src_cen(3));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(1)-src_cen(1));
            b_arr=[bij-dx bij bij+dx bij];
 
        elseif (plane_src==3)
            
            z=abs(obs_cen(1)-src_cen(1));
            aij=abs(obs_cen(2)-src_cen(2));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(3)-src_cen(3));
            b_arr=[bij-dx bij bij+dx bij];
            
        end
        a2minusz2over2=0.5*(a_arr.^2-z^2);
        b2minz2over2=0.5*(b_arr.^2-z^2);
        b2minus2z2=b_arr.^2-2*z^2;
        bz=b_arr.*z;
        
        int_res=0;
        for kk=1:4
            for mm=1:4
                rkm=sqrt(a_arr(kk)^2+b_arr(mm)^2+z^2);
                mult_fact=(-1)^(mm+kk);
                if abs(a_arr(kk)+rkm)<eps2
                    term1=0;
                else
                    
                    term1=b2minz2over2(mm)*a_arr(kk)*log(a_arr(kk)+rkm);
                end
                if abs(b_arr(mm)+rkm)<eps2
                    term2=0;
                else
                    
                    term2=a2minusz2over2(kk)*b_arr(mm)*log(b_arr(mm)+rkm);
                end
                
                term3=(1/6)*(b2minus2z2(mm)+a_arr(kk)^2)*rkm;
                if abs(z)<eps2
                    term4=0;
                else
                    
                    term4=bz(mm)*a_arr(kk)*atan((a_arr(kk)*b_arr(mm))/(rkm*z));
                end
                int_res=int_res+mult_fact*(term1+term2-term3-term4);
            end
        end
        
    elseif(fl_parallel == 0) % use the formula for othogonal panels
        if (plane_src == 1 && plane_obs == 2 || plane_src == 2 && plane_obs == 1)% xy,xz
            aij=abs(obs_cen(1)-src_cen(1));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(2)-src_cen(2));
            b_arr=[bij+dx/2 bij-dx/2];
            cij=abs(obs_cen(3)-src_cen(3));
            c_arr=[cij+dx/2 cij-dx/2];

        elseif (plane_src == 1 && plane_obs == 3 || plane_src == 3 && plane_obs == 1) % xy,yz
            aij=abs(obs_cen(2)-src_cen(2));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(1)-src_cen(1));
            b_arr=[bij+dx/2 bij-dx/2];
            cij=abs(obs_cen(3)-src_cen(3));
            c_arr=[cij+dx/2 cij-dx/2];

        elseif (plane_src == 2 && plane_obs == 3 || plane_src == 3 && plane_obs == 2)% xz,yz
            aij=abs(obs_cen(3)-src_cen(3));
            a_arr=[aij-dx aij aij+dx aij];
            bij=abs(obs_cen(1)-src_cen(1));
            b_arr=[bij+dx/2 bij-dx/2];
            cij=abs(obs_cen(2)-src_cen(2));
            c_arr=[cij+dx/2 cij-dx/2];

        end
        
        a2over2=a_arr.^2/2;
        c2over6=c_arr.^2/6;
        b2over6=b_arr.^2/6;
        a3over6=a_arr.^3/6;
        b2over2=b_arr.^2/2;
        c3over2=c_arr.^2/2;
        
        int_res=0;
        for kk=1:4
            for mm=1:2
                for ll=1:2
                    rkml=sqrt(a_arr(kk)^2+b_arr(mm)^2+c_arr(ll)^2);
                    mult_fact=(-1)^(mm+kk+ll+1);
                    if abs(b_arr(mm)+rkml)<eps2
                        term1=0;
                    else

                        term1=(a2over2(kk)-c2over6(ll))*c_arr(ll)*log(b_arr(mm)+rkml);
                    end
                    if abs(c_arr(ll)+rkml)<eps2
                        term2=0;
                    else
                 
                        term2=(a2over2(kk)-b2over6(mm))*b_arr(mm)*log(c_arr(ll)+rkml);
                    end
                    if abs(a_arr(kk)+rkml)<eps2
                        term3=0;
                    else
                        term3=a_arr(kk)*b_arr(mm)*c_arr(ll)*log(a_arr(kk)+rkml);
                    end
                    term4=b_arr(mm)*c_arr(ll)*rkml/3;
                    % in case of sigularity, use if-case
                    if abs(rkml*a_arr(kk))<eps2
                        term5=0;
                    else
                        
                        term5=a3over6(kk)*atan(b_arr(mm)*c_arr(ll)/(rkml*a_arr(kk)));
                    end
                    if abs(rkml*b_arr(mm))<eps2
                        term6=0;
                    else
                        
                        term6=b2over2(mm)*a_arr(kk)*atan((a_arr(kk)*c_arr(ll))/(rkml*b_arr(mm)));
                    end
                    if abs(c_arr(ll)*rkml)<eps2
                        term7=0;
                    else
                        
                        term7=a_arr(kk)*c3over2(ll)*atan((a_arr(kk)*b_arr(mm))/(c_arr(ll)*rkml));
                    end
                    int_res=int_res+mult_fact*(term1+term2+term3-term4-term5-term6-term7);
                end
            end
        end
    end
end





