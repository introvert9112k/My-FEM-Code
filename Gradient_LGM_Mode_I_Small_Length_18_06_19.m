% Using Small Large Length Scale Factor %
% Using Normalized Anisotropic Gradient Matrix %

% Gradient_LGM_Mode_I_Small_Length_18_06_19

%-----Localizing Gradient Damage Model based on Smooothed Stresses and
%Anisotropic Interaction Kernel------% 
clc;
tic;
close all;
L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 80; % Number of elements in X direction
numy = 80; % Number of elements in Y direction
stressState='PLANE_STRAIN'; %This defines the stressState chosen

%---------------------Material Parameters------------------%
nu = 0.2;  % Poisson's Ratio
E  = 1000; % Elastic Moduli,in MPa
cc = 0.04*L; % Internal length i.e. Small Length Scale paramter L/25
h = 1;
material_p =[nu, E, cc, h]; % Material Parameters

% Newton-Raphson Parameters

nsteps = 30; % Number of Load Increments i.e. Load is applied in "nsteps" increments
tol = 0.00001; % Tolerance Value to check Newton-Raphson convergence
maxit = 100; % Maximum Number of iterations in a single Load increment step i.e. Newton-Raphson iterations

%-----------------Damage Parameters---------------------%
K = 10; % Sensitivity Parameter i.e. Governs the sensitivity in compression relative to that in tension. 
       % Usually it is set equal to the ratio of compressive strength and tensile strength.
       %For calculating the modified von moises equivalent strain.
alpha = 0.99; % Controls the residual interactions,present in damage evolution law
beta = 9; % Controls the Slope of Softening Curve,present in damage evolution law.
hcoup = E*1e-9; % Coupling Modulus
R = 0.04; %require in interaction function g.
eta = 4; %require in interaction function g.
ft = 2;  % Tensile Strength in MPa
damage_p =[K alpha beta hcoup R eta]; % Damage Parameters


%-------------------------------Stiffness Matrix D--------------------------------------%
%If the stae of stress is plane stress
%De is the stiffness matrix for the plane stress case.
if ( strcmp(stressState,'PLANE_STRESS') )
    De = E/(1-nu^2)*[ 1   nu 0 ;
                      nu  1  0 ;
                       0   0  0.5*(1-nu) ];
%If the state of stress is plane strain,
%De is the stiffness matrix for the plane strain case.
elseif ( strcmp(stressState,'PLANE_STRAIN') )
    De = E/((1+nu)*(1-2*nu))*[1-nu  nu  0;
                              nu  1-nu 0;
                              0     0  (1-2*nu)/2];
end

%----------------------Element Type Selection and Mesh Generation-----------------------%

disp([num2str(toc),'   MESH GENERATION'])

%coordinates of the four corners of the plate
pt1 = [0 0] ; pt2 = [L 0] ; pt3 = [L D] ; pt4 = [0 D] ;
elemType1 = 'Q4' ;
elemType2 = 'Q4' ;

%node1 stores the location of all the nodes in global coordinate system.
%element1 stores the row numbers corresponding to the nodes of particular element in the node1. 
%For example node numbers for the 1st element is [1,2,83,82].That mean
%coordinates at these locations in the node1 matrix form element1.
[node1,element1] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType1);

[node2,element2] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType2);
numnode2 = size(node2,1);

if ( strcmp(elemType1,'Q8') )  
    % ------------------For Q8 Elements-------------------%
    % Four corner points
    nnx = 2*numx+1; % Number of nodes in X-direction
    nny = 2*numy+1; % Number of nodes in Y-direction
    numnode1 = size(node1,1);
    numelem = size(element1,1);
    urn =(nnx*nny)-(numx*numy);% upper right node number
    uln =urn-(nnx-1); % upper left node number
    lrn = nnx; % lower right node number
    lln = 1; % lower left node number
    
    % GET NODES ON ESSENTIAL BOUNDARY
    topEdge = [ uln:1:(urn-1); (uln+1):1:urn ]';
    botEdge = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
    botNodes = unique(botEdge);
    topNodes = unique(topEdge);
    dispNodes = botNodes; % Displacement B.C nodes
    tracNodes = topNodes; % Traction B.C. nodes
        
    % GET NODES For MODE I Problem
    rnodes = 1:(2*(numx/2));
    dispNodes1 = setdiff(dispNodes,rnodes);

elseif ( strcmp(elemType1,'Q4') ) 
    % ------------------For Q4 Elements-------------------%
    nnx = numx+1; % Number of nodes in X-direction
    nny = numy+1; % Number of nodes in Y-direction
    numnode1 = size(node1,1);
    numelem = size(element1,1);
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
    
    % GET NODES For MODE I Problem
    rnodes = 1:(numx/2);
    dispNodes1 = setdiff(dispNodes,rnodes);
end

%-----------------------Plot Mesh-------------------------------% 
disp([num2str(toc),'  PLOTTING MESH'])
figure
hold on
cntr = plot([0,L,L,0,0],[0,0,D,D,0],'k-');
set(cntr,'LineWidth',2);
plot_mesh(node1,element1,elemType1,'k-');
set(gcf, 'color', 'white');
plot(node1(dispNodes,1),node1(dispNodes,2),'ks');
plot(node1(tracNodes,1),node1(tracNodes,2),'ks');
plot(node1(:,1),node1(:,2),'.');
axis equal
axis off


%----------------Fixed DOF's for Boundary Conditions-----------------------%
udofs = 2*tracNodes(1,1) - 1;
vdofs = dispNodes1.*2;

upper_disp_node = tracNodes.*2;
lower_nodes = (dispNodes)*2;
l_nodes = [dispNodes.*2; dispNodes.*2-1];

%--------------------------------------------------------------------------%

%-----------------------Input Parameters-----------------------------% 
total_disp = numnode1*2; % Total unknown displacements % 2 dof per node
total_strain = numnode2*3; % Total unknown strains i.e. 3 additional dof per node
total_unknown = total_disp + total_strain; % Total Unknowns
u_tot =  zeros(total_disp,1); % Setting total unknown displacements to 0 initially
strain_tot = zeros(total_strain,1); % Setting total nonlocal strains to 0 initially
%numelem mean no of elements in the mesh
kappa = zeros(4*numelem,1); % History Parameter, value of kappa at each gauss point.
D_st= zeros(4*numelem,1); % Damage Variable vector,damage value at each gauss point.
NE_gp = zeros(4*numelem,1); % % Non-equivalent strain at each gauss point
stress_gp = zeros(4*numelem,3); %stress values at each gauss point.We have 3 stresses at each gauss point.
forcevdisp = zeros(2,nsteps+1);
forcevdisp(1,1) = 0;
forcevdisp(2,1) = 0;

[Gpnt] = computing_gauss_location(numelem,elemType2,node2,element2); % Gauss Point Location Generation

kappa0 = zeros(4*numelem,1); % Setting Initial History Variable
kappa0(:,1) = ft/E; % Tensile Strength/Elastic Modulus 'ft/E'


DAMAGE_DATA = []; %damage data
NESTRAIN_DATA = []; %non local equivalent strain data
DISP_DATA = []; %displacement data
NESTRAIN_DATA_NODES = []; %non local equivalent strain at nodes
INTERNAL_FORCE = []; %internal force data
INTERACTION_DATA = [];
SIGMA_XX = []; 
SIGMA_YY = [];
SIGMA_XY = [];    
SIGMA_XX_smooth = [];
SIGMA_YY_smooth = [];
SIGMA_XY_smooth = [];
EQ_STRESS = []; %equivalent stress
NEQ_STRESS = []; %non equivalent stress

%-----------------------Newton Raphson Loop-----------------------------% 
disp([num2str(toc),'  NEWTON RAPHSON LOOP BEGINS'])
for step = 1 : nsteps
    err3 = 1; %this is for intial tolerance.
    nit = 0;  %current iterations for this step
    Fint = zeros(total_unknown,1);
    fprintf(1,'\n Step %f \n',step);
    
    %-----------For 270 Load Steps-----------%   
    if step<=5 
    ubar = 0.012; %displacement for the steps <= 5    
    elseif step>5 
    ubar = 3e-3; %displacement for the steps > 5.
    end
    
    %Iterate until either the answer is converged or max iterations are 
    %reached.
    while ((err3>tol) && (nit<maxit))          % Newton Raphson loop
       
        nit = nit + 1; %incrementing the iterations
        
        % Computing Stiffness Matrix and Internal Force Vector
        [K,FAI,FE,D_st,kappa,NE_gp,stress_gp,interaction,stress_gp_sm,eq_stress,neq_stress] = globalstiffness(u_tot,strain_tot,D_st,material_p,De,damage_p,numelem,total_disp,...
            total_strain,node1,element1,element2,node2,elemType1,elemType2,kappa0,kappa,NE_gp,stress_gp);                     
             
        Fint = [FAI;FE]; %internal forces vector
        R = Fint;
        
        % ------------------------------------------------
        %       Imposing Essential Boundary Conditions
        %-------------------------------------------------
%         disp([num2str(toc),'   IMPOSING ESSENTIAL BOUNDARY CONDITION'])        
        
        Kres = K;
        bcwt = 1; % Used to keep the conditioning of the K matrix
        Kres(udofs,:) = 0;   % zero out the rows and columns of the K matrix
        Kres(:,udofs) = 0;
        Kres(vdofs,:) = 0;
        Kres(:,vdofs) = 0;
        Kres(udofs,udofs) = bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
        Kres(vdofs,vdofs) = bcwt*speye(length(vdofs));        
        Kres(upper_disp_node,:) = 0;
        Kres(:,upper_disp_node) = 0;
        Kres(upper_disp_node,upper_disp_node) = bcwt*speye(length(upper_disp_node));
        
        if nit == 1 
            for kk = 1:total_unknown               
                for ll = 1:length(upper_disp_node)
                    zz = upper_disp_node(ll);
                    R(kk,:) = R(kk,:) - K(kk,zz)*ubar;
                end
            end
        end
        
        R(udofs) = 0 ; % Imposing B.C. on Residual
        R(vdofs) = 0 ; % Imposing B.C. on Residual
        
        if nit == 1
            R(upper_disp_node) = ubar;
        else
            R(upper_disp_node) = 0;
        end
        
        %------Solve for the correction---------%
        [du]= Kres\R;
        du1 = du(1:total_disp,1); % Displacement increment vector
        du2 = du((total_disp+1):total_unknown,1); % Non-Equivalent Strain increment vector
        u_tot = u_tot + du1; % Updating displacement vector
        strain_tot = strain_tot + du2; % Updating strain vector                
        %------Checking convergence-----------------%

        wnorm1 = dot(u_tot,u_tot);
        wnorm2 = dot(strain_tot,strain_tot);
        err1 = dot(du1,du1);
        err2 = dot(du2,du2);
        err3 = dot(R,R);
        err1 = sqrt(err1/wnorm1); % rms error
        err2 = sqrt(err2/wnorm2);
        err3 = sqrt((err3)/(2*numnode1));
        fprintf(1,'Iteration number %d Correction-u %f Correction-ne %f Residual %f tolerance %f\n',nit,err1,err2,err3,tol);
    end
       
        DAMAGE_DATA(:,step) = D_st;
        NESTRAIN_DATA(:,step) = NE_gp;
        GPT_DATA = Gpnt;
        DISP_DATA(:,step) = u_tot;
        NESTRAIN_DATA_NODES(:,step) = strain_tot;
        INTERNAL_FORCE(:,step) = Fint;
        SIGMA_XX(:,step) = stress_gp(:,1);
        SIGMA_YY(:,step) = stress_gp(:,2);
        SIGMA_XY(:,step) = stress_gp(:,3);       
        SIGMA_XX_smooth(:,step) = stress_gp_sm(:,1);
        SIGMA_YY_smooth(:,step) = stress_gp_sm(:,2);
        SIGMA_XY_smooth(:,step) = stress_gp_sm(:,3);
        INTERACTION_DATA(:,step) = interaction;
        EQ_STRESS(:,step) = eq_stress;
        NEQ_STRESS(:,step) = neq_stress;
        
        forcevdisp(1,step+1) = mean(u_tot(upper_disp_node,:));
        forcevdisp(2,step+1) = sum(Fint(lower_nodes,:));   
        forcevdisp(3,step+1) = sum(Fint(l_nodes,:)); 
        
        save('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9.mat','DAMAGE_DATA','NESTRAIN_DATA','GPT_DATA','forcevdisp','DISP_DATA','NESTRAIN_DATA_NODES','INTERNAL_FORCE',...
            'INTERACTION_DATA','SIGMA_XX','SIGMA_YY','SIGMA_XY','SIGMA_XX_smooth','SIGMA_YY_smooth','SIGMA_XY_smooth','EQ_STRESS','NEQ_STRESS');
end

disp([num2str(toc),'  END OF NEWTON RAPHSON LOOP'])

% End of MaiN PrograM COde functioN
