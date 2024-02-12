function plot_local_equivalent_stress
load('Mode_I_steps_20_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_cycle.mat');
% subplot dimension

steps = [1,2,3,4,5; 6,7,8,9,10; 11,12,13,14,15; 16,17,18,19,20];
n1 = size(steps,1); % number of rows
n2 = size(steps,2); % number of columns

guass_point_locations = GPT_DATA;
local_equivalent_stress = EQ_STRESS;

figure;
for k1 = 1:n1
    for k2 = 1:n2
            st = steps(k1,k2);
            subplot(n1,n2,(k1-1)*n2 + k2);
            tri = delaunay(guass_point_locations(:,1),guass_point_locations(:,2));
            patch('Vertices',guass_point_locations,'Faces',tri,'FaceVertexCData',local_equivalent_stress(:,st));  
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
sgtitle('Local Equivalent stress');
end 