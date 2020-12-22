function [precond,ids]=implement_preconditioner_ver2(dx,num_vox_in_blk,L,M,N,geom_bndry_panels_test,Eps_inout,inds_glob,num_diel)
lili=tic;
fl_diel=0;
if (sum(abs(Eps_inout(:,1))) > 1e-13) % check whether we have a dielectric panel
    fl_diel = 1; % yes, we have
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fl_diel==1)
    
    x_aligned_cond_inds=find(abs(geom_bndry_panels_test(:,4)) == 1 & ...
        geom_bndry_panels_test(:,7) == 0);
    x_aligned_diel_inds=find(abs(geom_bndry_panels_test(:,4)) == 1 & ...
        geom_bndry_panels_test(:,7) > 0);
    
    y_aligned_cond_inds=find(abs(geom_bndry_panels_test(:,4)) == 2 & ...
        geom_bndry_panels_test(:,7) == 0);
    y_aligned_diel_inds=find(abs(geom_bndry_panels_test(:,4)) == 2 & ...
        geom_bndry_panels_test(:,7) > 0);
    
    z_aligned_cond_inds=find(abs(geom_bndry_panels_test(:,4)) == 3 & ...
        geom_bndry_panels_test(:,7) == 0);
    z_aligned_diel_inds=find(abs(geom_bndry_panels_test(:,4)) == 3 & ...
        geom_bndry_panels_test(:,7) > 0);
    
    inds_sorted=[x_aligned_cond_inds; x_aligned_diel_inds; ...
        y_aligned_cond_inds; y_aligned_diel_inds; z_aligned_cond_inds; z_aligned_diel_inds];
    
    geom_bndry_panels_test=geom_bndry_panels_test(inds_sorted,:);
    clear x_aligned_cond_inds x_aligned_diel_inds y_aligned_cond_inds y_aligned_diel_inds z_aligned_cond_inds z_aligned_diel_inds inds_sorted
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num_big_blk_x=ceil(L/num_vox_in_blk);  num_big_blk_y=ceil(M/num_vox_in_blk);  num_big_blk_z=ceil(N/num_vox_in_blk);
% num=num_vox_in_blk;

num_blk_x=ceil(L/num_vox_in_blk);  num_blk_y=ceil(M/num_vox_in_blk);  num_blk_z=ceil(N/num_vox_in_blk);
if mod(L,num_vox_in_blk)>0 && L>1
    num_big_blk_x=num_blk_x-1;
else
    num_big_blk_x=num_blk_x;
end
if mod(M,num_vox_in_blk)>0 && M>1
    num_big_blk_y=num_blk_y-1;
else
    num_big_blk_y=num_blk_y;
end
if mod(N,num_vox_in_blk)>0 && N>1
    num_big_blk_z=num_blk_z-1;
else
    num_big_blk_z=num_blk_z;
end

num=max([L-(floor(L/num_vox_in_blk)-1)*num_vox_in_blk,M-(floor(M/num_vox_in_blk)-1)*num_vox_in_blk,N-(floor(N/num_vox_in_blk)-1)*num_vox_in_blk]);



xinxin=tic;
% ids=cell(num_big_blk_x*num_big_blk_y*num_big_blk_z,1);
% % parpool('local',31);
% parfor tt=1:num_big_blk_x*num_big_blk_y*num_big_blk_z
%     ids{tt}=find(geom_bndry_panels_test(:,8)==tt);
% end
% ids(cellfun(@isempty,ids))=[];

ids=cell(num_big_blk_x*num_big_blk_y*num_big_blk_z,1);
% a=geom_bndry_panels_test(:,8);
[b,c]=sort(geom_bndry_panels_test(:,10));
num_ind=zeros(num_big_blk_x*num_big_blk_y*num_big_blk_z,1);
for kk=1:length(geom_bndry_panels_test(:,10))
    num_ind(b(kk))=num_ind(b(kk))+1;
end
clear b;
st_ind=1;
for ll=1:num_big_blk_x*num_big_blk_y*num_big_blk_z
    end_ind=st_ind+num_ind(ll)-1;
    ids{ll}=c(st_ind:end_ind);
    st_ind=end_ind+1;
end
clear c;
ids(cellfun(@isempty,ids))=[];
xin=toc(xinxin);
disp(['precond_new_ids ::: ',num2str(xin)])
% ids=ids1;
% ids(cellfun(@isempty,ids))=[];
%%% -----------------------------------------------------------------------------------------------------
%%%                              Basic elements in preconditioner
%%% ------------------------------------------------------------------------------------------------------
% inds_glob(1,:)=[1 num_x_aligned_cond];
% inds_glob(num_diel+2,:) = [num_x_aligned_all+1 num_x_aligned_all+num_y_aligned_cond];
% inds_glob(num_diel*2+3,:) = [num_x_aligned_all+num_y_aligned_all+1 num_x_aligned_all+num_y_aligned_all+num_z_aligned_cond]; % conductor part

geom_bndry_panels_test(inds_glob(num_diel+2,1):inds_glob(num_diel+2,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel+2,1):inds_glob(num_diel+2,2),1:3)+[-dx/2,dx/2,0];
geom_bndry_panels_test(inds_glob(num_diel*2+3,1):inds_glob(num_diel*2+3,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel*2+3,1):inds_glob(num_diel*2+3,2),1:3)+[-dx/2,0,dx/2];
if num_diel>0
    for ii=1:num_diel
        geom_bndry_panels_test(inds_glob(num_diel+2+ii,1):inds_glob(num_diel+2+ii,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel+2+ii,1):inds_glob(num_diel+2+ii,2),1:3)+[-dx/2,dx/2,0];
        geom_bndry_panels_test(inds_glob(num_diel*2+3+ii,1):inds_glob(num_diel*2+3+ii,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel*2+3+ii,1):inds_glob(num_diel*2+3+ii,2),1:3)+[-dx/2,0,dx/2];
    end
end
% geom_bndry_panels_test(find(abs(geom_bndry_panels_test(:,4))==2),1:3)=geom_bndry_panels_test(find(abs(geom_bndry_panels_test(:,4))==2),1:3)+[-dx/2,dx/2,0];
% geom_bndry_panels_test(find(abs(geom_bndry_panels_test(:,4))==3),1:3)=geom_bndry_panels_test(find(abs(geom_bndry_panels_test(:,4))==3),1:3)+[-dx/2,0,dx/2];

tic
[fN_charge]=precond_circulant(dx,num,num,num,Eps_inout(end,:));
disp(['precond_circulate ::: ',num2str(toc)])

precond_test=cell(size(ids,1),1);
titi=tic;
% parpool('local',31);
geom_bndry_panels_test(:,8:10)=[];
geom_bndry_panels_test(:,5:6)=[];
parfor ii=1:size(ids,1)
    precond_test{ii}=zeros(length(ids{ii}),length(ids{ii}));
%     precond_test{ii}=precond_retrieval(fN_charge,ids{ii},geom_bndry_panels_test,dx,num);
    index=zeros(1,3);
    precond_test{ii}=zeros(length(ids{ii}),length(ids{ii}));
    for ll=1:length(ids{ii})
        for kk=1:length(ids{ii})
            if kk==ll
                if geom_bndry_panels_test(ids{ii}(kk),5)==0
                    precond_test{ii}(kk,ll)=fN_charge{1,1,1}(1,1,1);
                else
                    precond_test{ii}(kk,ll)=fN_charge{1,1,2}(1,1,1);
                end
%             elseif kk < ll
%                 if geom_bndry_panels_test(ids{ii}(kk),5)==0
%                     if abs(geom_bndry_panels_test(ids{ii}(kk),4))==1 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==1
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{1,1,1}(index(1),index(2),index(3));
%                     elseif abs(geom_bndry_panels_test(ids{ii}(kk),4))==2 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==2
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{2,2,1}(index(1),index(2),index(3));
%                     elseif abs(geom_bndry_panels_test(ids{ii}(kk),4))==3 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==3
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{3,3,1}(index(1),index(2),index(3));
%                     elseif abs(geom_bndry_panels_test(ids{ii}(kk),4))==1 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==2
%                         
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{1,2,1}(index(1),index(2),index(3));
%                     elseif abs(geom_bndry_panels_test(ids{ii}(kk),4))==1 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==3
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{1,3,1}(index(1),index(2),index(3));
%                     elseif abs(geom_bndry_panels_test(ids{ii}(kk),4))==2 && abs(geom_bndry_panels_test(ids{ii}(ll),4))==3
%                         if geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1)>=-1e-13
%                             index(1)=round((geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1);
%                         else
%                             index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids{ii}(kk),1)-geom_bndry_panels_test(ids{ii}(ll),1))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2)>=-1e-13
%                             index(2)=round((geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1);
%                         else
%                             index(2)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),2)-geom_bndry_panels_test(ids{ii}(ll),2))/dx+1)+1);
%                         end
%                         if geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3)>=-1e-13
%                             index(3)=round((geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1);
%                         else
%                             index(3)=round((2*num+(geom_bndry_panels_test(ids{ii}(kk),3)-geom_bndry_panels_test(ids{ii}(ll),3))/dx+1)+1);
%                         end
%                         precond_test{ii}(kk,ll)=fN_charge{2,3,1}(index(1),index(2),index(3));
%                     end
%                 end
            end
        end
    end
end

% parfor ii=1:size(ids,1);precond_test{ii}=precond_test{ii}+(precond_test{ii})'+diag(diag(precond_test{ii}));end %-diag(diag(precond_test{ii}))
toto=toc(titi);
disp(['precond_retrieval ::: ',num2str(toto)])


% ind=zeros(size(ids,1),1);
% for ii=1:size(ids,1);ind(ii,1)=length(ids{ii});end
% [~,I1]=sort(ind);
% ids=ids(I1);
% precond_test=precond_test(I1);



tic
% LL=cell(size(ids,1),1);
% UU=cell(size(ids,1),1);
% PP=cell(size(ids,1),1);
precond=cell(size(ids,1),1);
% da=cell(size(ids,1),1);
% parpool('local',31);
parfor kk=1:size(ids,1)
    if sum(ids{kk})>0
%         [LL,UU,PP]=lu(precond_test{kk});
        precond{kk}=inv(precond_test{kk});
        precond_test{kk}=[];
%         precond{kk}=UU\(LL\(PP*eye(size(ids{kk},1))));
        %         precond{kk}=inv(precond_test{kk});
        %         da=decomposition(precond_test{kk})
        %         precond{kk}=da\eye(size(ids{kk},1));
    end
end
disp(['precond_inverse ::: ',num2str(toc)])

time_precond=toc(lili);
disp(['Total time for generating preconditioner ::: ' ,num2str(time_precond)]);

average_unk_in_blk=size(geom_bndry_panels_test,1)/size(ids,1);
disp(['The average unknowns in each preconditioner block ::: ' ,num2str(average_unk_in_blk)]);
