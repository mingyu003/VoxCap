function [JOut_full]=matvect_mult_charge_pecdiel_multi_diel1(ch_coef_dum, fN_ch_all,inds_glob,locs_diel_pos,ids_panel_diel,ids_panels,L,M,N,num,ss,fl_Tucker_decomp,ids,precond,seg,ord,fl_precond,ind_layer)
% the function is used to mutiply circulent tensor with charge density to
% generate potential vector.
%--------------------------------------------------------------------------
%%%INPUTS:
%ch_coef_dum: charge density;
%fN_compression: circulant tensor after compression and the first cell is the size of these tensors;
%inds_glob: intervals from start panel to end panel in x, y and z
%           direction.
%locs_diel_pos: locations of elements of diel_pnl_ids
%ids_panels: the numbering of conductor panels in 3 directions.
%---------------------------------------------------------------------------
% 0) put the charge coefficients in seperate vectors due to their panel
% alignments
num_cond=num(1);num_diel=num(2);num_interface=num(3);
tot_num=sum(num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      right precond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if fl_precond==2
%     for ll=1:size(ids,1)
%         if size(precond{ll},1)==2 && size(precond{ll},2)==1
%             ch_coef_dum(ids{ll},:)=precond{ll}(1)*ch_coef_dum(ids{ll},:);
%         else
%             ch_coef_dum(ids{ll},:)=precond{ll}*ch_coef_dum(ids{ll},:);
%         end
%     end
% elseif fl_precond==1
%     for kk=1:length(I)
%         for ll=1:size(seg,1)
%             if kk>=seg(ll,1)-1e-3 && kk<=seg(ll,2)+1e-3
%                 if size(precond{ll},1)==2 && size(precond{ll},2)==1
%                     ch_coef_dum(ids{kk},:)=precond{ll}(1)*ch_coef_dum(ids{kk},:);
%                 else
%                     ch_coef_dum(ids{I(kk)},:)=precond{ll}*ch_coef_dum(ids{I(kk)},:);
%                 end
%             end
%         end
%     end
% end
if fl_precond==2
    for ll=1:size(ids,1)
        ch_coef_dum(ids{ll},:)=precond{ll}*ch_coef_dum(ids{ll},:);
    end
    if isempty(ord)==0
        for kk=1:size(ord,1)
            for ii=1:size(ord{kk})
%                 if ss(ind_layer{kk}(ii))>1e-13
                    ch_coef_dum(ord{kk}(ii),:)=(1/ss(ind_layer{kk}(ii)))*ch_coef_dum(ord{kk}(ii),:);
%                 end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (num_diel > 0)
    ch_coefs_xx = ch_coef_dum(inds_glob(1,1):inds_glob(tot_num,2)); % [rhoa_x_c;rhoa_x_d]
    ch_coefs_yy = ch_coef_dum(inds_glob(tot_num+1,1):inds_glob(2*(tot_num),2)); % [rhoa_y_c;rhoa_y_d]
    ch_coefs_zz = ch_coef_dum(inds_glob(2*(tot_num)+1,1):inds_glob(3*(tot_num),2)); % [rhoa_z_c;rhoa_z_d]
    
    coe_xx=cell(num_diel+num_interface,1);coe_yy=cell(num_diel+num_interface,1);coe_zz=cell(num_diel+num_interface,1);
    
    for ii=1:num_diel
        coe_xx{ii}=ss(ii)*ch_coef_dum(inds_glob(num_cond+ii,1):inds_glob(num_cond+ii,2)); % coefficient for self-term extraction
        coe_yy{ii}=ss(ii)*ch_coef_dum(inds_glob(num_cond+tot_num+ii,1):inds_glob(num_cond+tot_num+ii,2));
        coe_zz{ii}=ss(ii)*ch_coef_dum(inds_glob(num_cond+2*tot_num+ii,1):inds_glob(num_cond+2*tot_num+ii,2));
    end
    for ii=1:num_interface
%         coe_xx{num_diel+ii}=ss(num_diel+ii)*ch_coef_dum(inds_glob(num_cond+num_diel+ii,1):inds_glob(num_cond+num_diel+ii,2)); % coefficient for self-term extraction
%         coe_yy{num_diel+ii}=ss(num_diel+ii)*ch_coef_dum(inds_glob(num_cond+num_diel+tot_num+ii,1):inds_glob(num_cond+num_diel+tot_num+ii,2));
        coe_zz{num_diel+ii}=ss(num_diel+ii)*ch_coef_dum(inds_glob(num_cond+num_diel+2*tot_num+ii,1):inds_glob(num_cond+num_diel+2*tot_num+ii,2));
    end
    
%     coe_xx1=ss(1)*ch_coef_dum(inds_glob(2,1):inds_glob(2,2)); % coefficient for self-term extraction
%     coe_yy1=ss(1)*ch_coef_dum(inds_glob(5,1):inds_glob(5,2));
%     coe_zz1=ss(1)*ch_coef_dum(inds_glob(8,1):inds_glob(8,2));
%     
%     coe_xx2=ss(2)*ch_coef_dum(inds_glob(3,1):inds_glob(3,2)); % coefficient for self-term extraction
%     coe_yy2=ss(2)*ch_coef_dum(inds_glob(6,1):inds_glob(6,2));
%     coe_zz2=ss(2)*ch_coef_dum(inds_glob(9,1):inds_glob(9,2));
else
    ch_coefs_xx = ch_coef_dum(inds_glob(1,1):inds_glob(num_cond,2)); % [rhoa_x_c]
    ch_coefs_yy = ch_coef_dum(inds_glob(num_cond+1,1):inds_glob(2*num_cond,2)); % [rhoa_y_c]
    ch_coefs_zz = ch_coef_dum(inds_glob(2*num_cond+1,1):inds_glob(3*num_cond,2)); % [rhoa_z_c]
end

JOut_full=zeros(size(ch_coef_dum,1),1);

% 1) get the sizes of FFTs
if fl_Tucker_decomp==1
LfN_x = fN_ch_all{1}(1,1); MfN_x = fN_ch_all{1}(1,2); NfN_x = fN_ch_all{1}(1,3);  %xx
LfN_xy = fN_ch_all{1}(4,1); MfN_xy = fN_ch_all{1}(4,2); NfN_xy = fN_ch_all{1}(4,3); %xy
LfN_xz = fN_ch_all{1}(5,1); MfN_xz = fN_ch_all{1}(5,2); NfN_xz = fN_ch_all{1}(5,3); %xz
else
[LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
[LfN_xy, MfN_xy, NfN_xy, ~] = size(fN_ch_all{1,2,1}); %xy
[LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
end
% 2) put the charge coefficients in tensors
JIn_x = zeros(L+1, M, N);
temp=[];
for ii=1:tot_num-num_interface
% for ii=1:tot_num-1
    temp=[temp;ids_panels{ii}{1}];
end
JIn_x(temp) = ch_coefs_xx;

JIn_y = zeros(L, M+1, N);
temp=[];
% for ii=1:tot_num-1
for ii=1:tot_num-num_interface
    temp=[temp;ids_panels{ii}{2}];
end
JIn_y(temp) = ch_coefs_yy;

JIn_z = zeros(L, M, N+1);
temp=[];
for ii=1:tot_num
    temp=[temp;ids_panels{ii}{3}];
end
JIn_z(temp) = ch_coefs_zz;

% fftw('dwisdom',[]); 
% fftw('planner','measure');

fJ_x = fftn(JIn_x(:,:,:),[LfN_x, MfN_x, NfN_x]); % fJ_yx = fJ_zx =fJ_x;
clear JIn_x;
fJ_y = fftn(JIn_y(:,:,:),[LfN_xy, MfN_xy, NfN_xy]); % fJ_zy = fJ_xy = fJ_y;
clear JIn_y;
fJ_z = fftn(JIn_z(:,:,:),[LfN_xz, MfN_xz, NfN_xz]); % fJ_yz = fJ_xz = fJ_z;
clear JIn_z;

% conductor
if fl_Tucker_decomp==1
fN_ch_c_xy=ten_mat_prod(fN_ch_all{14},{fN_ch_all{15},fN_ch_all{16},fN_ch_all{17}});%ten_mat_prod
Jout1_x=fJ_y .* fN_ch_c_xy;
Jout1_y=fJ_x .* conj(fN_ch_c_xy);
clear fN_ch_c_xy;
fN_ch_c_xz=ten_mat_prod(fN_ch_all{18},{fN_ch_all{19},fN_ch_all{20},fN_ch_all{21}});
Jout1_x=Jout1_x+fJ_z .* fN_ch_c_xz;
Jout1_z=fJ_x .* conj(fN_ch_c_xz);
clear fN_ch_c_xz;
fN_ch_c_yz=ten_mat_prod(fN_ch_all{22},{fN_ch_all{23},fN_ch_all{24},fN_ch_all{25}});
Jout1_y=Jout1_y+fJ_z .* fN_ch_c_yz;
Jout1_z=Jout1_z+fJ_y .* conj(fN_ch_c_yz);
clear fN_ch_c_yz;
Jout1_x = Jout1_x + fJ_x .* ten_mat_prod(fN_ch_all{2},{fN_ch_all{3},fN_ch_all{4},fN_ch_all{5}});
Jout1_x=ifftn(Jout1_x);

Jout1_y=Jout1_y+fJ_y .* ten_mat_prod(fN_ch_all{6},{fN_ch_all{7},fN_ch_all{8},fN_ch_all{9}});
Jout1_y=ifftn(Jout1_y);

Jout1_z=Jout1_z+fJ_z .* ten_mat_prod(fN_ch_all{10},{fN_ch_all{11},fN_ch_all{12},fN_ch_all{13}});
Jout1_z=ifftn(Jout1_z);


else

    Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,1}+fJ_y .* fN_ch_all{1,2,1}+fJ_z .* fN_ch_all{1,3,1});% xx,xy,xz
    Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,1}+fJ_x .* conj(fN_ch_all{1,2,1})+fJ_z .* fN_ch_all{2,3,1});% yx,yy,yz
    Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,1}+fJ_x .* conj(fN_ch_all{1,3,1})+fJ_y .* conj(fN_ch_all{2,3,1}));% zx,zy,zz
   
end    

Jout1_x = Jout1_x(1:L+1,1:M,1:N);
for ii=1:tot_num-num_interface
    JOut_full(inds_glob(ii,1):inds_glob(ii,2)) = Jout1_x(ids_panels{ii}{1});
end

Jout1_y = Jout1_y(1:L,1:M+1,1:N);
for ii=tot_num+1:2*tot_num-num_interface
    JOut_full(inds_glob(ii,1):inds_glob(ii,2)) = Jout1_y(ids_panels{ii-tot_num}{2});
end

Jout1_z = Jout1_z(1:L,1:M,1:N+1);
for ii=2*tot_num+1:3*tot_num
    if isempty(ids_panels{ii-2*tot_num}{3})==0
    JOut_full(inds_glob(ii,1):inds_glob(ii,2)) = Jout1_z(ids_panels{ii-2*tot_num}{3});
    end
end

if (num_diel > 0)
    % dielectric
    if fl_Tucker_decomp==1
    Jout1_x = ifftn(fJ_x .* ten_mat_prod(fN_ch_all{38},{fN_ch_all{39},fN_ch_all{40},fN_ch_all{41}}) + ...
        fJ_y .* ten_mat_prod(fN_ch_all{50},{fN_ch_all{51},fN_ch_all{52},fN_ch_all{53}}) + ...
        fJ_z .* ten_mat_prod(fN_ch_all{54},{fN_ch_all{55},fN_ch_all{56},fN_ch_all{57}})); % xx,xy,xz
    Jout1_y = ifftn(fJ_y .* ten_mat_prod(fN_ch_all{42},{fN_ch_all{43},fN_ch_all{44},fN_ch_all{45}})+ ...
        fJ_x .* ten_mat_prod(fN_ch_all{62},{fN_ch_all{63},fN_ch_all{64},fN_ch_all{65}})+ ...
        fJ_z .* ten_mat_prod(fN_ch_all{58},{fN_ch_all{59},fN_ch_all{60},fN_ch_all{61}})); %yx,yy,yz
    Jout1_z = ifftn(fJ_z .* ten_mat_prod(fN_ch_all{46},{fN_ch_all{47},fN_ch_all{48},fN_ch_all{49}})+ ...
        fJ_x .* ten_mat_prod(fN_ch_all{66},{fN_ch_all{67},fN_ch_all{68},fN_ch_all{69}})+ ...
        fJ_y .* ten_mat_prod(fN_ch_all{70},{fN_ch_all{71},fN_ch_all{72},fN_ch_all{73}}));  % zx,zy,zz
    else

       Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,2} + fJ_y .* fN_ch_all{1,2,2} + fJ_z .* fN_ch_all{1,3,2});% xx,xy,xz
       Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,2}+fJ_x .* fN_ch_all{2,1,2}+fJ_z .* fN_ch_all{2,3,2});%yx,yy,yz
       Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,2}+fJ_x .* fN_ch_all{3,1,2}+fJ_y .* fN_ch_all{3,2,2});   % zx,zy,zz
    end
    
       
    Jout1_x = Jout1_x(1:L+1,1:M,1:N);
    for ii=1:num_diel
        JOut_full(inds_glob(num_cond+ii,1)+locs_diel_pos{ii}{1}-1) = Jout1_x(ids_panel_diel{ii}{1});
        JOut_full(inds_glob(num_cond+ii,1)+locs_diel_pos{ii}{4}-1) = -Jout1_x(ids_panel_diel{ii}{4});
        JOut_full(inds_glob(num_cond+ii,1):inds_glob(num_cond+ii,2))=coe_xx{ii}+JOut_full(inds_glob(num_cond+ii,1):inds_glob(num_cond+ii,2));
    end
%     for ii=1:num_interface
%         JOut_full(inds_glob(num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{1}-1) = Jout1_x(ids_panel_diel{ii+num_diel}{1});
%         JOut_full(inds_glob(num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{4}-1) = -Jout1_x(ids_panel_diel{ii+num_diel}{4});
%         JOut_full(inds_glob(num_cond+num_diel+ii,1):inds_glob(num_cond+num_diel+ii,2))=coe_xx{num_diel+ii}+JOut_full(inds_glob(num_cond+num_diel+ii,1):inds_glob(num_cond+num_diel+ii,2));
%     end
    
    Jout1_y = Jout1_y(1:L,1:M+1,1:N);
    for ii=1:num_diel
        JOut_full(inds_glob(tot_num+num_cond+ii,1)+locs_diel_pos{ii}{2}-1) = Jout1_y(ids_panel_diel{ii}{2});
        JOut_full(inds_glob(tot_num+num_cond+ii,1)+locs_diel_pos{ii}{5}-1) = -Jout1_y(ids_panel_diel{ii}{5});
        JOut_full(inds_glob(tot_num+num_cond+ii,1):inds_glob(tot_num+num_cond+ii,2))=coe_yy{ii}+JOut_full(inds_glob(tot_num+num_cond+ii,1):inds_glob(tot_num+num_cond+ii,2));
    end
%     for ii=1:num_interface
%         JOut_full(inds_glob(tot_num+num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{2}-1) = Jout1_y(ids_panel_diel{ii+num_diel}{2});
%         JOut_full(inds_glob(tot_num+num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{5}-1) = -Jout1_y(ids_panel_diel{ii+num_diel}{5});
%         JOut_full(inds_glob(tot_num+num_cond+num_diel+ii,1):inds_glob(tot_num+num_cond+num_diel+ii,2))=coe_yy{ii+num_diel}+JOut_full(inds_glob(tot_num+num_cond+num_diel+ii,1):inds_glob(tot_num+num_cond+num_diel+ii,2));
%     end
   
    Jout1_z = Jout1_z(1:L,1:M,1:N+1);
    for ii=1:num_diel
        JOut_full(inds_glob(2*tot_num+num_cond+ii,1)+locs_diel_pos{ii}{3}-1) = Jout1_z(ids_panel_diel{ii}{3});
        JOut_full(inds_glob(2*tot_num+num_cond+ii,1)+locs_diel_pos{ii}{6}-1) = -Jout1_z(ids_panel_diel{ii}{6});
        JOut_full(inds_glob(2*tot_num+num_cond+ii,1):inds_glob(2*tot_num+num_cond+ii,2))=coe_zz{ii}+JOut_full(inds_glob(2*tot_num+num_cond+ii,1):inds_glob(2*tot_num+num_cond+ii,2));
    end
    for ii=1:num_interface
        JOut_full(inds_glob(2*tot_num+num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{3}-1) = Jout1_z(ids_panel_diel{ii+num_diel}{3});
        JOut_full(inds_glob(2*tot_num+num_cond+num_diel+ii,1)+locs_diel_pos{ii+num_diel}{6}-1) = -Jout1_z(ids_panel_diel{ii+num_diel}{6});
        JOut_full(inds_glob(2*tot_num+num_cond+num_diel+ii,1):inds_glob(2*tot_num+num_cond+num_diel+ii,2))=coe_zz{ii+num_diel}+JOut_full(inds_glob(2*tot_num+num_cond+num_diel+ii,1):inds_glob(2*tot_num+num_cond+num_diel+ii,2));
    end
end

% fprintf ('.') ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      left precond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if fl_precond==1
%     for ll=1:size(ids,1)
%         JOut_full(ids{ll},:)=precond{ll}*JOut_full(ids{ll},:);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
