function precond_test=precond_retrieval(fN_charge,ids,geom_bndry_panels_test,dx,num)
index=zeros(1,3);
precond_test=zeros(length(ids),length(ids));
for ll=1:length(ids)
    for kk=1:length(ids)
        if kk==ll
            if geom_bndry_panels_test(ids(kk),5)==0
                precond_test(kk,ll)=fN_charge{1,1,1}(1,1,1);
            else
                precond_test(kk,ll)=fN_charge{1,1,2}(1,1,1);
            end
        else%if kk < ll
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