%--------------------------Plot Output:Damage-----------------------------%
function plot_damage_alok

check_step = 30;
%--------------------- For Single Damage Plot-------------------------%
% figure
% hold on
% tri = delaunay(GPT_DATA(:,1),GPT_DATA(:,2));
% patch('Vertices',GPT_DATA,'Faces',tri,'FaceVertexCData',DAMAGE_DATA(:,check_step));
% % plot_mesh(node,element,elemType,'r-');
% % plot(GPT_DATA(:,1),GPT_DATA(:,2),'.');
% colormap('jet');
% colorbar 
% shading interp
% set(gcf, 'color', 'white');
% axis equal
% axis off

load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Tension_cycle.mat');
% subplot dimension
n1 = 2; % number of rows
steps = [1,2,3,4,5; 6,7,8,9,10];
n2 = 5; % number of columns

guass_point_locations = GPT_DATA;
damage = DAMAGE_DATA;

for k1 = 1:n1
    for k2 = 1:n2
            st = steps(k1,k2);
            % subplot(n1,n2,(k1-1)*n2 + k2,...
            % 'position', [(1-nw)/n2/2 + (k2-1)/n2, (1-nh)/n1/2 + 1-k1/n1,...
            % nw/n2 nh/n1]);
            subplot(n1,n2,(k1-1)*n2 + k2);
            tri = delaunay(guass_point_locations(:,1),guass_point_locations(:,2));
            patch('Vertices',guass_point_locations,'Faces',tri,'FaceVertexCData',damage(:,st));  
            %plot_mesh(node,element,elemType,'r-');
            hold on
            colormap('jet');
            shading interp
            set(gcf, 'color', 'white');
            axis equal
            axis off
            title(sprintf("load Step %d"), (k1-1)*n2 + k2);
            colorbar;
    end
end
end 

% End of the function 
