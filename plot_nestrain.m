%--------------------------Plot Output:Damage-----------------------------%
%-----------------Plot Output:Non-Equivalent Strain---------------------------%

function plot_nestrain

close all; clc;
L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 80; % Number of elements in X direction
numy = 80; % Number of elements in Y direction
%load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9'); 
load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Compression.mat');

% Principle Stress Based Localizing GDM
% step = [80 165 270]; % Small Length Scale
step = [10];  % Large Length Scale

% Conventional Localizing GDM
% step = [30 75 230]; % Small Length Scale
% step = [22 75 200];  % Large Length Scale

loading = 'MODE_I';

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
% u_x = DISP_DATA(1:2:2*numnode,step1);
% u_y = DISP_DATA(2:2:2*numnode,step1);
% updated_node(:,1) = node(:,1) + u_x;
% updated_node(:,2) = node(:,2) + u_y;

%--------------------- Damage Plot-------------------------%
figure

% subplot dimension
n1 = 1; % number of rows
n2 = 3; % number of columns

% These values would define the space between the graphs
% if equal to 1 there will be no space between graphs
nw = 0.95; % normalized width
nh = 0.95; % normalized height

for k1 = 1:n1
    for k2 = 1:n2
            check_step = step(k2);
            subplot(n1,n2,(k1-1)*n2 + k2,...
            'position', [(1-nw)/n2/2 + (k2-1)/n2, (1-nh)/n1/2 + 1-k1/n1,...
            nw/n2 nh/n1]);
            tri = delaunay(GPT_DATA(:,1),GPT_DATA(:,2));
            patch('Vertices',GPT_DATA,'Faces',tri,'FaceVertexCData',NESTRAIN_DATA(:,check_step));
            hold on
            plot([0.5 L/2],[0.5 0.5],'-k','LineWidth',2);            
%             plot_mesh(node,element,elemType,'r-');
            colormap('jet');
            shading interp
            set(gcf, 'color', 'white');
            axis equal
            axis off
    end
end

end

% End of the function