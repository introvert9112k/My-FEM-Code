function plot_non_local_eq_strain_gp


% loading = 'MODE_I';
% % loading = 'COMPRESSION';
% 
% %------------------Material Parameters------------------%
% if (strcmp(loading,'MODE_I') )
%     E  = 1000; % Elastic Moduli
%     kappa0 = 0.002;
% elseif (strcmp(loading,'COMPRESSION') )
%     E  = 20000; % Elastic Moduli
%     kappa0 = 0.0001;
% end



%----------------------------------At Single load Step---------------------------
% check_step = 30;
% figure
% hold on
% tri = delaunay(GPT_DATA(:,1),GPT_DATA(:,2));
% patch('Vertices',GPT_DATA,'Faces',tri,'FaceVertexCData',NESTRAIN_DATA(:,check_step));
% % plot_mesh(node,element,elemType,'r-');
% colormap('jet');
% colorbar 
% shading interp
% set(gcf, 'color', 'white');
% axis equal
% axis off
% % title("Non Local Equivalent Strain");
% end 


%--------------------- For multiple load steps-----------------------
%--------------------- Micro-morphic-equivalent-strain-plot--------------------------%
figure;
%load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Tension_cycle.mat');
%load('Mode_I_steps_20_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_cycle.mat');
load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_cycle_0.012.mat');
% subplot dimension

%steps = [1,2,3,4,5; 6,7,8,9,10; 11,12,13,14,15; 16,17,18,19,20];
steps = [1,2,3,4,5; 6,7,8,9,10];
n1 = size(steps,1); % number of rows
n2 = size(steps,2); % number of columns

guass_point_locations = GPT_DATA;
micromorphic_equivalent_strain = NESTRAIN_DATA;

for k1 = 1:n1
    for k2 = 1:n2
            st = steps(k1,k2);
            % subplot(n1,n2,(k1-1)*n2 + k2,...
            % 'position', [(1-nw)/n2/2 + (k2-1)/n2, (1-nh)/n1/2 + 1-k1/n1,...
            % nw/n2 nh/n1]);
            subplot(n1,n2,(k1-1)*n2 + k2);
            tri = delaunay(guass_point_locations(:,1),guass_point_locations(:,2));
            patch('Vertices',guass_point_locations,'Faces',tri,'FaceVertexCData',micromorphic_equivalent_strain(:,st));  
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
 sgtitle('Non local Equivalent Strain');
end  
