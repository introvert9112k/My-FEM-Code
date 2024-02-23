function plot_force_vs_displacement
close all; clc;
L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 60; % Number of elements in X direction
numy = 60; % Number of elements in Y direction

loading = 'Tension';
% loading = 'Compression';
%load('Mode_I_steps_20_1_by_1_Eta_4_R04_SmallLenScale_Beta_9.mat');
load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9.mat');

%------------------Material Parameters------------------%
if (strcmp(loading,'Tension') )
    E  = 1000; % Elastic Moduli
    kappa0 = 0.002;
elseif (strcmp(loading,'Compression') )
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
%---------------------Load Displacement Plot---------------------%
figure
hold on
%plot((forcevdisp(1,:)/L)*1e3,forcevdisp(2,:)/(L*E*kappa0),'-r','LineWidth',1);
plot( forcevdisp(1,:),forcevdisp(2,:),'-r','LineWidth',1);
xlabel({'(u/L)*1e3'},'FontSize',12);
ylabel({'F/(L*E*K0)'},'FontSize',12);

end 