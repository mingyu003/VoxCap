function [precond,ids,seg]=implement_preconditioner(dx,num_vox_in_blk,L,M,N,geom_bndry_panels_test,Eps_inout,inds_glob,num_diel)
% This subroutine genrates preconditioner with almost all unique block. the
% memory requirement is relatively small, but the method to judge whether
% two blocks are same is not mathematical. After some testings,if the judge
% method is satisfied, the two blocks are the same but if  the judge method is 
% not satisfied, the two blocks may be the same. 

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

num_big_blk_x=ceil(L/num_vox_in_blk);  num_big_blk_y=ceil(M/num_vox_in_blk);  num_big_blk_z=ceil(N/num_vox_in_blk);
num=num_vox_in_blk;

tic
[fN_charge]=precond_circulant(dx,num,num,num,Eps_inout(end,:));
disp(['precond_circulate ::: ',num2str(toc)])


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

size_ids=size(ids,1);
ind=zeros(size_ids,1);
for ii=1:size(ids,1);ind(ii,1)=length(ids{ii});end
[B,I]=sort(ind);
ids=ids(I);

geom_bndry_panels_test(:,8:10)=[];
geom_bndry_panels_test(:,5:6)=[];
geom_bndry_panels_test(inds_glob(num_diel+2,1):inds_glob(num_diel+2,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel+2,1):inds_glob(num_diel+2,2),1:3)+[-dx/2,dx/2,0];
geom_bndry_panels_test(inds_glob(num_diel*2+3,1):inds_glob(num_diel*2+3,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel*2+3,1):inds_glob(num_diel*2+3,2),1:3)+[-dx/2,0,dx/2];
if num_diel>0
    for ii=1:num_diel
        geom_bndry_panels_test(inds_glob(num_diel+2+ii,1):inds_glob(num_diel+2+ii,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel+2+ii,1):inds_glob(num_diel+2+ii,2),1:3)+[-dx/2,dx/2,0];
        geom_bndry_panels_test(inds_glob(num_diel*2+3+ii,1):inds_glob(num_diel*2+3+ii,2),1:3)=geom_bndry_panels_test(inds_glob(num_diel*2+3+ii,1):inds_glob(num_diel*2+3+ii,2),1:3)+[-dx/2,0,dx/2];
    end
end


prec=cell(3,1);
for ii=1:4
    prec{ii}=precond_retrieval(fN_charge,ids{ii},geom_bndry_panels_test,dx,num);
end

if size(prec{2},1)==size(prec{1},1) 
    if norm(prec{2}-prec{1})/norm(prec{2})<1e-13
        aa=ids{2}-ids{1};
    elseif size(prec{3},1)==size(prec{1},1)
        if  norm(prec{3}-prec{1})/norm(prec{3})<1e-13
            aa=ids{3}-ids{1};
        end
    end
elseif size(prec{2},1)~=size(prec{1},1) && size(prec{3},1)==size(prec{2},1) 
    if norm(prec{2}-prec{3})/norm(prec{2})<1e-13 
    aa=ids{3}-ids{2};
    end
end

clear prec;

factors=[aa(1) aa(end)]; 
clear aa;

diff_num=zeros(length(unique(B)),1);
count=1; diff_num(1)=1;
for ii=1:size_ids-1

    if (abs(B(ii+1)-B(ii))>1e-13)
        count=count+1;
        diff_num(count)=ii+1;
    end
end

counter=0;
sameid=cell(size(ids,1),1);
segment=zeros(size(ids,1),2);
for ii=1:count-1
    diff_vect=diff_num(ii):diff_num(ii+1)-1;
    for nn=1:diff_num(ii+1)-diff_num(ii)
%         st_num=1;
        id_same=[diff_vect(1)];
        
        for mm=2:length(diff_vect)
            ll=diff_vect(mm);
            aa=ids{ll}-ids{diff_vect(1)};
            if mod(aa(1),factors(1))==0 && aa(end)>0 && mod(aa(end),factors(2))==0
                id_same=[id_same;ll];
            end
        end
        counter=counter+1;
        sameid{counter}=id_same;
        diff_vect=setdiff(diff_vect,id_same);
        segment(counter,1:2)=[1,length(id_same)];
        if length(diff_vect)==0
            break;
        end
    end
end
% ii=count
diff_vect=diff_num(count):size_ids;
if length(diff_vect)>1
    for nn=1:length(diff_vect)
        %         st_num=1;
        id_same=[diff_vect(1)];
        
        for mm=2:length(diff_vect)
            ll=diff_vect(mm);
            aa=ids{ll}-ids{diff_vect(1)};
            if mod(aa(1),factors(1))==0 && aa(end)>0 && mod(aa(end),factors(2))==0
                id_same=[id_same;ll];
            end
        end
        counter=counter+1;
        sameid{counter}=id_same;
        diff_vect=setdiff(diff_vect,id_same);
        segment(counter,1:2)=[1,length(id_same)];
        if length(diff_vect)==0
            break;
        end
    end
else
    sameid{counter+1}=[diff_num(end):size(ids,1)]';
    segment(counter+1,1:2)=[1 length(sameid{counter+1})];
end
%%%%%%%%%%%%%%%
sameid(cellfun(@isempty,sameid))=[];  
segment(all(segment==0,2),:)=[]; 
inds_ids=[];
for ii=1:size(sameid,1)
    inds_ids=[inds_ids;sameid{ii}];
end
ids=ids(inds_ids);

seg=zeros(size(sameid,1),2);
seg(1,1:2)=segment(1,1:2);
for ii=2:size(sameid,1)
    seg(ii,1:2)=[seg(ii-1,2)+1,seg(ii-1,2)+1+segment(ii,2)-segment(ii,1)];
end

          
clear B

%%% -----------------------------------------------------------------------------------------------------
%%%                              Basic elements in preconditioner
%%% ------------------------------------------------------------------------------------------------------

count=size(sameid,1);

precond_test=cell(count,1);

% c=parcluster;
% p=c.parpool(12);

titi=tic;
% parpool('local',31);
% I=1:size(ids,1);
parfor jj=1:count
    ii=seg(jj,1);   
    precond_test{jj}=precond_retrieval(fN_charge,ids{ii},geom_bndry_panels_test,dx,num);
end

% parfor ii=1:count
%     if size(precond_test{ii},1)~=2 && size(precond_test{ii},2)~=1
%     precond_test{ii}=precond_test{ii}+(precond_test{ii})';
%     end
% end
toto=toc(titi);
disp(['precond_retrieval ::: ',num2str(toto)])

tic
% LL=cell(size(ids,1),1);
% UU=cell(size(ids,1),1);
% PP=cell(size(ids,1),1);
precond=cell(count,1);
% da=cell(size(ids,1),1);
% parpool('local',31);
parfor kk=1:count
    if sum(ids{kk})>0
        %         [LL,UU,PP]=lu(precond_test{kk});
        if size(precond_test{kk},1)==2 && size(precond_test{kk},2)==1
            precond{kk}=[1/precond_test{kk}(1);precond_test{kk}(2)];
        else
            precond{kk}=inv(precond_test{kk});
            precond_test{kk}=[];
        end
        %         precond{kk}=UU\(LL\(PP*eye(size(ids{kk},1))));
        %         precond{kk}=inv(precond_test{kk});
        %         da=decomposition(precond_test{kk})
        %         precond{kk}=da\eye(size(ids{kk},1));
    end
end
disp(['precond_inverse ::: ',num2str(toc)])
% p.delete;
time_precond=toc(lili);
disp(['Total time for generating preconditioner ::: ' ,num2str(time_precond)]);

average_unk_in_blk=size(geom_bndry_panels_test,1)/size(ids,1);
disp(['The average unknowns in each preconditioner block ::: ' ,num2str(average_unk_in_blk)]);
