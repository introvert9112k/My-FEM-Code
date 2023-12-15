function plot_displacement ()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load('Mode_I_steps_20_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_cycle.mat');
% subplot dimension
n1 = 4; % number of rows
steps = [1,2,3,4,5; 6,7,8,9,10; 11,12,13,14,15; 16,17,18,19,20];
n2 = 5; % number of columns

displacement_data = DISP_DATA;
nRows = size(displacement_data,1);
nodes_locations = node1;
vertical_displacement_data = zeros(nRows/2,n1*n2);

for k1 = 1:n1
    for k2 = 1:n2
        st = steps(k1,k2);
        for i = 2 : 2 :nRows 
            vertical_displacement_data(i/2,st) = DISP_DATA(i,st);
        end 
    end
end 


for k1 = 1:n1
    for k2 = 1:n2
            st = steps(k1,k2);
            subplot(n1,n2,(k1-1)*n2 + k2);
            tri = delaunay(nodes_locations(:,1),nodes_locations(:,2));
            patch('Vertices',nodes_locations,'Faces',tri,'FaceVertexCData',vertical_displacement_data(:,st));  
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