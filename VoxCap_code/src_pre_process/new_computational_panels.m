function [unique_panel]=new_computational_panels(dx,grid_intcon)
% the function is used to generate computational domain
%fl_check_pnls=1, plot all panels



[L,M,N,~]=size(grid_intcon);
%%%------------------------------------------------------------------
%             obtain total panels' coordinates and normal
%%%------------------------------------------------------------------
total_panels=zeros(L,M,N,24);
for kk=1:L
    for ll=1:M
        for mm=1:N
            total_panels(kk,ll,mm,1)=grid_intcon(kk,ll,mm,1)-dx/2;
            total_panels(kk,ll,mm,2)=grid_intcon(kk,ll,mm,2);
            total_panels(kk,ll,mm,3)=grid_intcon(kk,ll,mm,3);
            total_panels(kk,ll,mm,4)=1;
            total_panels(kk,ll,mm,9)=grid_intcon(kk,ll,mm,1);
            total_panels(kk,ll,mm,10)=grid_intcon(kk,ll,mm,2)-dx/2;
            total_panels(kk,ll,mm,11)=grid_intcon(kk,ll,mm,3);
            total_panels(kk,ll,mm,12)=2;
            total_panels(kk,ll,mm,17)=grid_intcon(kk,ll,mm,1);
            total_panels(kk,ll,mm,18)=grid_intcon(kk,ll,mm,2);
            total_panels(kk,ll,mm,19)=grid_intcon(kk,ll,mm,3)-dx/2;
            total_panels(kk,ll,mm,20)=3;
            if kk==L
            total_panels(kk,ll,mm,5)=grid_intcon(kk,ll,mm,1)+dx/2;
            total_panels(kk,ll,mm,6)=grid_intcon(kk,ll,mm,2);
            total_panels(kk,ll,mm,7)=grid_intcon(kk,ll,mm,3);
            total_panels(kk,ll,mm,8)=1;
            end
            if ll==M
            total_panels(kk,ll,mm,13)=grid_intcon(kk,ll,mm,1);
            total_panels(kk,ll,mm,14)=grid_intcon(kk,ll,mm,2)+dx/2;
            total_panels(kk,ll,mm,15)=grid_intcon(kk,ll,mm,3);
            total_panels(kk,ll,mm,16)=2;  
            end
            if mm==N
            total_panels(kk,ll,mm,21)=grid_intcon(kk,ll,mm,1);
            total_panels(kk,ll,mm,22)=grid_intcon(kk,ll,mm,2);
            total_panels(kk,ll,mm,23)=grid_intcon(kk,ll,mm,3)+dx/2;
            total_panels(kk,ll,mm,24)=3; 
            end
        end
    end
end
total_panels=reshape(total_panels,L*M*N,24);
unique_panel=zeros(L*M*N*6,4);
unique_panel(1:L*M*N,:)=total_panels(:,1:4);
unique_panel((1+L*M*N):2*L*M*N,:)=total_panels(:,5:8);
unique_panel((1+2*L*M*N):3*L*M*N,:)=total_panels(:,9:12);
unique_panel((1+3*L*M*N):4*L*M*N,:)=total_panels(:,13:16);
unique_panel((1+4*L*M*N):5*L*M*N,:)=total_panels(:,17:20);
unique_panel((1+5*L*M*N):6*L*M*N,:)=total_panels(:,21:24); 
unique_panel(all(unique_panel==0,2),:)=[];

% unique_panel=unique(unique_panel,'rows'); 
% clear total_panels;
%%%-------------------------------------------------------------------
%                  obtain numbering
%%%-------------------------------------------------------------------
idsx=1:(L+1)*M*N;
tempx=idsx((L+1):(L+1):end);
idsx((L+1):(L+1):end)=0;
idsx(idsx==0)=[];
idsx=[idsx' ; tempx'];
clear tempx;

idsy=1:(M+1)*L*N;
tempy=cell(1,20);
for mm=1:N
    tempy{mm}=idsy((L*M+L*(M+1)*(mm-1)+1):(L*M+L*(M+1)*(mm-1)+L));
    idsy((L*M+L*(M+1)*(mm-1)+1):(L*M+L*(M+1)*(mm-1)+L))=0;
end
idsy(idsy==0)=[];
tempy=cell2mat(tempy);
idsy=[idsy' ; tempy'];
idsy=idsy+(L+1)*M*N;
clear tempy;

idsz=1:L*M*(N+1);
idsz=idsz';
idsz=idsz+(M+1)*L*N+(L+1)*M*N;

ids=[idsx;idsy;idsz];
unique_panel(:,5)=ids;
clear idsx idsy idsz;
unique_panel(:,6)=1;
% %%%----------------------------------------------------------------
% %                        plot structure
% %%%----------------------------------------------------------------
% if (fl_check_pnls > 0)
%     
%     % select what you'll plot - all panels or only boundary panels - one of
%     % the following:
%     if (fl_check_pnls  == 1)
%         tmp_mat=unique_panel;
%     end
%     
%     figure
%     for kk=1:size(tmp_mat,1)
%         verts=[];
%         if (tmp_mat(kk,4) == 1) % x-directed
%             verts=[tmp_mat(kk,1) tmp_mat(kk,2)-dx/2 tmp_mat(kk,3)-dx/2; ...
%                 tmp_mat(kk,1) tmp_mat(kk,2)+dx/2 tmp_mat(kk,3)-dx/2; ...
%                 tmp_mat(kk,1) tmp_mat(kk,2)+dx/2 tmp_mat(kk,3)+dx/2; ...
%                 tmp_mat(kk,1) tmp_mat(kk,2)-dx/2 tmp_mat(kk,3)+dx/2;];
%         elseif (tmp_mat(kk,4) == 2) % y-directed
%             verts=[tmp_mat(kk,1)-dx/2 tmp_mat(kk,2) tmp_mat(kk,3)-dx/2; ...
%                 tmp_mat(kk,1)+dx/2 tmp_mat(kk,2) tmp_mat(kk,3)-dx/2; ...
%                 tmp_mat(kk,1)+dx/2 tmp_mat(kk,2) tmp_mat(kk,3)+dx/2; ...
%                 tmp_mat(kk,1)-dx/2 tmp_mat(kk,2) tmp_mat(kk,3)+dx/2;];
%         elseif (tmp_mat(kk,4) == 3) % z-directed
%             verts=[tmp_mat(kk,1)-dx/2 tmp_mat(kk,2)-dx/2 tmp_mat(kk,3); ...
%                 tmp_mat(kk,1)+dx/2 tmp_mat(kk,2)-dx/2 tmp_mat(kk,3); ...
%                 tmp_mat(kk,1)+dx/2 tmp_mat(kk,2)+dx/2 tmp_mat(kk,3); ...
%                 tmp_mat(kk,1)-dx/2 tmp_mat(kk,2)+dx/2 tmp_mat(kk,3);];
%         end
%         patch(verts(:,1),verts(:,2),verts(:,3),[1 2 4 3],'EdgeColor','blue','FaceColor','none');
%         hold on
%         %plot3(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),'r+')
%         %hold on
%         %text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(kk,5))
%         text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(tmp_mat(kk,5)))
%         %text(tmp_mat(kk,1),tmp_mat(kk,2),tmp_mat(kk,3),num2str(tmp_mat(kk,6)))
%         hold on
%     end
%     xlabel('x');ylabel('y');zlabel('z');view(-30,30);set(gca,'FontSize',16); axis tight;
% end


disp('Done... Extracting unique panels of a grid')
disp('-----------------------------------------------------')

