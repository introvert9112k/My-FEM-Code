function plot_non_local_strain_nodes
% -----------------------Description------------------------------
% Contour Plot of the Vertical Micromorphic strain at each Node
load('Mode_I_steps_20_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_cycle.mat');
%load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Tension_cycle.mat');
% subplot dimension

steps = [1,2,3,4,5; 6,7,8,9,10; 11,12,13,14,15;16,17,18,19,20];
%steps = [1,2,3,4,5; 6,7,8,9,10];
n1 = size(steps,1); % number of rows
n2 = size(steps,2); % number of columns

non_local_strains = NESTRAIN_DATA_NODES; 
nRows = size(non_local_strains,1);
node_locations = node1;
vertical_non_local_strain = zeros(nRows/3,n1*n2);

for k1 = 1:n1
    for k2 = 1:n2
        st = steps(k1,k2);
        ind = 1;
        for i = 2 : 3 : nRows 
            vertical_non_local_strain(ind,st) = non_local_strains(i,st);
            ind = ind + 1;
        end 
    end
end 

for k1 = 1:n1
    for k2 = 1:n2
            st = steps(k1,k2);
            subplot(n1,n2,(k1-1)*n2 + k2);
            tri = delaunay(node_locations(:,1),node_locations(:,2));
            patch('Vertices',node_locations,'Faces',tri,'FaceVertexCData',vertical_non_local_strain(:,st));  
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