%--------------------------Plot Output:Damage-----------------------------%

function plot_damage

close all; clc;
L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 60; % Number of elements in X direction
numy = 60; % Number of elements in Y direction
% load('Mode_I_80by80_Beta_60_30_10_18_3_eta_5_R_dot04');
load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9');

check_step = 1;

loading = 'MODE_I';
% loading = 'COMPRESSION';

%------------------Material Parameters------------------%
if (strcmp(loading,'MODE_I') )
    E  = 1000; % Elastic Moduli
    kappa0 = 0.002;
elseif (strcmp(loading,'COMPRESSION') )
    E  = 20000; % Elastic Moduli
    kappa0 = 0.0001;
end

pt1 = [0 0] ; pt2 = [L 0] ; pt3 = [L D] ; pt4 = [0 D] ;
elemType = 'Q4' ;


if ( strcmp(elemType,'Q8') )  
    [node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType);
    % ------------------For Q8 Elements-------------------%
    % Four corner points
    nnx = 2*numx+1; % Number of nodes in X-direction
    nny = 2*numy+1; % Number of nodes in Y-direction
    numnode = size(node,1);
    numelem = size(element,1);
    urn = (nnx*nny)-(numx*numy);% upper right node number
    uln = urn-(nnx-1); % upper left node number
    lrn = nnx; % lower right node number
    lln = 1; % lower left node number
    
    % GET NODES ON ESSENTIAL BOUNDARY
    topEdge = [ uln:1:(urn-1); (uln+1):1:urn ]';
    botEdge = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
    botNodes = unique(botEdge);
    topNodes = unique(topEdge);
    dispNodes = botNodes; % Displacement B.C nodes
    tracNodes = topNodes; % Traction B.C. nodes

elseif (strcmp(elemType,'Q4') ) 
    % ------------------For Q4 Elements-------------------%
    [node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType);
    nnx = numx+1; % Number of nodes in X-direction
    nny = numy+1; % Number of nodes in Y-direction
    numnode = size(node,1);
    numelem = size(element,1);
    uln = nnx*(nny-1)+1;       % upper left node number
    urn = nnx*nny;             % upper right node number
    lrn = nnx;                 % lower right node number
    lln = 1;                   % lower left node number
    
    % GET NODES ON ESSENTIAL BOUNDARY
    topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
    botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
    botNodes   = unique(botEdge);
    topNodes   = unique(topEdge);
    dispNodes = botNodes; % Displacement B.C nodes
    tracNodes = topNodes; % Traction B.C. nodes
end

% updated_node = zeros(numnode,2);
% u_x = DISP_DATA(1:2:2*numnode,check_step);
% u_y = DISP_DATA(2:2:2*numnode,check_step);
% updated_node(:,1) = node(:,1) + u_x;
% updated_node(:,2) = node(:,2) + u_y;

%--------------------- Damage Plot-------------------------%
figure
hold on
tri = delaunay(GPT_DATA(:,1),GPT_DATA(:,2));
patch('Vertices',GPT_DATA,'Faces',tri,'FaceVertexCData',DAMAGE_DATA(:,check_step));
% plot_mesh(node,element,elemType,'r-');
% plot(GPT_DATA(:,1),GPT_DATA(:,2),'.');
colormap('jet');
colorbar 
shading interp
set(gcf, 'color', 'white');
axis equal
axis off
title("Damage")

end 

% End of the function 
