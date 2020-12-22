function plot_panels(dx,plot_panels)
%%%----------------------------------------------------------------
%                        plot structure
%%%----------------------------------------------------------------


% select what you'll plot - all panels or only boundary panels - one of
% the following:

figure
for kk=1:size(plot_panels,1)
    verts=[];
    if (plot_panels(kk,4) == 1) % x-directed
        verts=[plot_panels(kk,1) plot_panels(kk,2)-dx/2 plot_panels(kk,3)-dx/2; ...
            plot_panels(kk,1) plot_panels(kk,2)+dx/2 plot_panels(kk,3)-dx/2; ...
            plot_panels(kk,1) plot_panels(kk,2)+dx/2 plot_panels(kk,3)+dx/2; ...
            plot_panels(kk,1) plot_panels(kk,2)-dx/2 plot_panels(kk,3)+dx/2;];
    elseif (plot_panels(kk,4) == 2) % y-directed
        verts=[plot_panels(kk,1)-dx/2 plot_panels(kk,2) plot_panels(kk,3)-dx/2; ...
            plot_panels(kk,1)+dx/2 plot_panels(kk,2) plot_panels(kk,3)-dx/2; ...
            plot_panels(kk,1)+dx/2 plot_panels(kk,2) plot_panels(kk,3)+dx/2; ...
            plot_panels(kk,1)-dx/2 plot_panels(kk,2) plot_panels(kk,3)+dx/2;];
    elseif (plot_panels(kk,4) == 3) % z-directed
        verts=[plot_panels(kk,1)-dx/2 plot_panels(kk,2)-dx/2 plot_panels(kk,3); ...
            plot_panels(kk,1)+dx/2 plot_panels(kk,2)-dx/2 plot_panels(kk,3); ...
            plot_panels(kk,1)+dx/2 plot_panels(kk,2)+dx/2 plot_panels(kk,3); ...
            plot_panels(kk,1)-dx/2 plot_panels(kk,2)+dx/2 plot_panels(kk,3);];
    end
    patch(verts(:,1),verts(:,2),verts(:,3),[1 2 4 3],'EdgeColor','blue','FaceColor','none');
    hold on
    text(plot_panels(kk,1),plot_panels(kk,2),plot_panels(kk,3),num2str(plot_panels(kk,5)))
    hold on
end
xlabel('x');ylabel('y');zlabel('z');view(-30,30);set(gca,'FontSize',16); axis tight;



disp('Done... Extracting unique panels of a grid')
disp('-----------------------------------------------------')