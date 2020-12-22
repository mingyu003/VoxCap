function [precond_test,ord,ind_layer]=precond_retrieval_multilayer(fN_charge,ids,geom_bndry_panels_test,dx,num)
index=zeros(1,3);
precond_test=zeros(length(ids),length(ids));
counter=1;
ord=zeros(length(ids),1);
ind_layer=zeros(length(ids),1);
for ll=1:length(ids)
    for kk=1:length(ids)
        if kk==ll
            if geom_bndry_panels_test(ids(kk),5)==0 %0: cond 1: diel
                precond_test(kk,ll)=fN_charge{1,1,1}(1,1,1);
            else
                if geom_bndry_panels_test(ids(kk),5)==1
                    ord(counter)=ids(kk);
                    ind_layer(counter)=1;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(1);
                elseif geom_bndry_panels_test(ids(kk),5)==2
                    ord(counter)=ids(kk);
                    ind_layer(counter)=2;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(2);
                elseif geom_bndry_panels_test(ids(kk),5)==3
                    ord(counter)=ids(kk);
                    ind_layer(counter)=3;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(3);
                elseif geom_bndry_panels_test(ids(kk),5)==4
                    ord(counter)=ids(kk);
                    ind_layer(counter)=4;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(4);
                elseif geom_bndry_panels_test(ids(kk),5)==5
                    ord(counter)=ids(kk);
                    ind_layer(counter)=5;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(5);
                elseif geom_bndry_panels_test(ids(kk),5)==6
                    ord(counter)=ids(kk);
                    ind_layer(counter)=6;
                    counter=counter+1;
%                     precond_test(kk,ll)=ss(6);
                elseif geom_bndry_panels_test(ids(kk),5)==7
                    ord(counter)=ids(kk);
                    ind_layer(counter)=7;
                    counter=counter+1;
                    %                     precond_test(kk,ll)=ss(7);
                elseif geom_bndry_panels_test(ids(kk),5)==8
                    ord(counter)=ids(kk);
                    ind_layer(counter)=8;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==9
                    ord(counter)=ids(kk);
                    ind_layer(counter)=9;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==10
                    ord(counter)=ids(kk);
                    ind_layer(counter)=10;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==11
                    ord(counter)=ids(kk);
                    ind_layer(counter)=11;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==12
                    ord(counter)=ids(kk);
                    ind_layer(counter)=12;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==13
                    ord(counter)=ids(kk);
                    ind_layer(counter)=13;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==14
                    ord(counter)=ids(kk);
                    ind_layer(counter)=14;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==15
                    ord(counter)=ids(kk);
                    ind_layer(counter)=15;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==16
                    ord(counter)=ids(kk);
                    ind_layer(counter)=16;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==17
                    ord(counter)=ids(kk);
                    ind_layer(counter)=17;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==18
                    ord(counter)=ids(kk);
                    ind_layer(counter)=18;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==19
                    ord(counter)=ids(kk);
                    ind_layer(counter)=19;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==20
                    ord(counter)=ids(kk);
                    ind_layer(counter)=20;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==21
                    ord(counter)=ids(kk);
                    ind_layer(counter)=21;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==22
                    ord(counter)=ids(kk);
                    ind_layer(counter)=22;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==23
                    ord(counter)=ids(kk);
                    ind_layer(counter)=23;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==24
                    ord(counter)=ids(kk);
                    ind_layer(counter)=24;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==25
                    ord(counter)=ids(kk);
                    ind_layer(counter)=25;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==26
                    ord(counter)=ids(kk);
                    ind_layer(counter)=26;
                    counter=counter+1;
                elseif geom_bndry_panels_test(ids(kk),5)==27
                    ord(counter)=ids(kk);
                    ind_layer(counter)=27;
                    counter=counter+1;
                end
            end
        elseif kk < ll
            if geom_bndry_panels_test(ids(kk),5)==0
                if abs(geom_bndry_panels_test(ids(kk),4))==1 && abs(geom_bndry_panels_test(ids(ll),4))==1
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*num+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{1,1,1}(index(1),index(2),index(3));
                elseif abs(geom_bndry_panels_test(ids(kk),4))==2 && abs(geom_bndry_panels_test(ids(ll),4))==2
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*num+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{2,2,1}(index(1),index(2),index(3));
                elseif abs(geom_bndry_panels_test(ids(kk),4))==3 && abs(geom_bndry_panels_test(ids(ll),4))==3
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*num+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{3,3,1}(index(1),index(2),index(3));
                elseif abs(geom_bndry_panels_test(ids(kk),4))==1 && abs(geom_bndry_panels_test(ids(ll),4))==2
                    
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*num+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*num+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{1,2,1}(index(1),index(2),index(3));
                elseif abs(geom_bndry_panels_test(ids(kk),4))==1 && abs(geom_bndry_panels_test(ids(ll),4))==3
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*num+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*num+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{1,3,1}(index(1),index(2),index(3));
                elseif abs(geom_bndry_panels_test(ids(kk),4))==2 && abs(geom_bndry_panels_test(ids(ll),4))==3
                    if geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1)>=-1e-13
                        index(1)=round((geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1);
                    else
                        index(1)=round((2*(num-1)+(geom_bndry_panels_test(ids(kk),1)-geom_bndry_panels_test(ids(ll),1))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2)>=-1e-13
                        index(2)=round((geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1);
                    else
                        index(2)=round((2*num+(geom_bndry_panels_test(ids(kk),2)-geom_bndry_panels_test(ids(ll),2))/dx+1)+1);
                    end
                    if geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3)>=-1e-13
                        index(3)=round((geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1);
                    else
                        index(3)=round((2*num+(geom_bndry_panels_test(ids(kk),3)-geom_bndry_panels_test(ids(ll),3))/dx+1)+1);
                    end
                    precond_test(kk,ll)=fN_charge{2,3,1}(index(1),index(2),index(3));
                end
            end
        end
    end
end

ord(counter:end,:)=[];
ind_layer(counter:end,:)=[];
% if sum(sum(precond_test))==0
%     dim=size(precond_test,1);
%     precond_test=[fN_charge{1,1,2}(1,1,1);dim];
% end

% for ll=1:length(ids)
%     for kk=1:length(ids)
%         if kk==ll
%             if geom_bndry_panels_test(ids(kk),5)==0
%                 precond_test(kk,ll)=fN_charge{1,1,1}(1,1,1);
%             else
%                 precond_test(kk,ll)=fN_charge{1,1,2}(1,1,1);
%             end
%         end
%     end
% end

