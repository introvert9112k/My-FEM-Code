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
nu = 0.3;  % Poisson's Ratio
E  = 210000; % Elastic Moduli in N/mm^2
cc = 0.01; % Internal length i.e. Small Length Scale,in mm
th = 0.5; % thickness of the plate, in mm
material_p = [nu, E, cc, th]; % Material Parameters

% Newton-Raphson Parameters
tol = 0.00001; % Tolerance Value to check Newton-Raphson convergence
maxit = 100; % Maximum Number of iterations in a single Load increment step i.e. Newton-Raphson iterations

%-----------------Damage Parameters---------------------%
K = 10; % Sensitivity Parameter i.e. Governs the sensitivity in compression relative to that in tension. 
       % Usually it is set equal to the ratio of compressive strength and tensile strength.
       %For calculating the modified von moises equivalent strain.
alpha = 10; % Controls the residual interactions
beta = 8.09; % Controls the Slope of Softening Curve
C = (6.60) * 10^21; 
hcoup = E*1e-9; % Coupling Modulus
R = 0.04; %require in interaction function g.
eta = 4; %require in interaction function g.
ft = 2;  % Tensile Strength in MPa
damage_p =[K alpha beta hcoup R eta,C]; % Damage Parameters

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
    dispNodes = botNodes; % Displacement B.C nodes [1,81]
    tracNodes = topNodes; % Traction B.C. nodes    [6481,6561]
    
    % GET NODES For MODE I Problem
    rnodes = 1:(numx/2);
    dispNodes1 = setdiff(dispNodes,rnodes); %Second half of the bottom edge nodes.
end

%-----------------------Plot Mesh-------------------------------% 
% disp([num2str(toc),'  PLOTTING MESH'])
% figure
% hold on
% cntr = plot([0,L,L,0,0],[0,0,D,D,0],'k-');
% set(cntr,'LineWidth',2);
% plot_mesh(node1,element1,elemType1,'k-');
% set(gcf, 'color', 'white');
% plot(node1(dispNodes,1),node1(dispNodes,2),'ks');
% plot(node1(tracNodes,1),node1(tracNodes,2),'ks');
% plot(node1(:,1),node1(:,2),'.');
% axis equal
% axis off


%----------------Fixed DOF's for Boundary Conditions-----------------------%
udofs = 2*tracNodes(1,1) - 1; %displacement along u-direction is not allowed for the first upper_node
vdofs = dispNodes1.*2; %displacement along the v-direction is not allowed for the lower_nodes

upper_disp_node = tracNodes.*2; %This gives the location of the vertical displacement for the topEdge Nodes in u_tot 
lower_nodes = (dispNodes)*2; %This gives the location of the vertical displacement for the bottomEdge Nodes in u_tot 
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


[Gpnt] = computing_gauss_location(numelem,elemType2,node2,element2); % Gauss Point Location Generation

kappa0 = zeros(4*numelem,1); % Setting Initial History Variable
kappa0(:,1) = 0.00114; % Tensile Strength/Elastic Modulus 'ft/E'


DAMAGE_DATA = []; %damage value at each guass point [25600] at each load step 30. size = [25600x30]  
NESTRAIN_DATA = []; %non local equivalent strain data at each guass point [25600] at each load step 30. size = [25600x30] 
DISP_DATA = []; %displacement at each node [81x81],each node has 2 displacements i.e 81x81x2 at each load step. size = [13122x30]
NESTRAIN_DATA_NODES = []; %Strain at nodes [81x81],each node has 3 strains i.e 81x81x3 at each load step. size = [19683x30]
INTERNAL_FORCE = []; %internal force data
INTERACTION_DATA = [];
SIGMA_XX = []; %sigma_xx stress at each guass point [25600] at each load step [25600x30]
SIGMA_YY = []; %sigma_xx stress at each guass point [25600] at each load step [25600x30]
SIGMA_XY = []; %sigma_xx stress at each guass point [25600] at each load step [25600x30]  
SIGMA_XX_smooth = []; 
SIGMA_YY_smooth = [];
SIGMA_XY_smooth = [];
EQ_STRESS = [];  %Equivalent stress at each guass point at each load step [25600x30]
NEQ_STRESS = []; %Non local Equivalent stress at each guass point at each load step [25600x30] 
STRAIN_LOCAL_XX = [];     %strain_xx stress at each guass point [25600] at each load step 
STRAIN_LOCAL_YY = [];     %strain_yy stress at each guass point [25600] at each load step 
STRAIN_LOCAL_XY  = [];    %strain_xy stress at each guass point [25600] at each load step 
STRAIN_NON_LOCAL_XX = []; %micormorphic or non local strain_xx at each guass point [25600] at each load step 
STRAIN_NON_LOCAL_YY = []; %micormorphic or non local strain_yy  stress at each guass point [25600] at each load step 
STRAIN_NON_LOCAL_XY = []; %micormorphic or non local strain_xy  at each guass point [25600] at each load step .
EQ_STRAIN = []; %Equivalent strain from local strain tensor.

%-----------------------Newton Raphson Loop-----------------------------% 
disp([num2str(toc),'  NEWTON RAPHSON LOOP BEGINS'])   


temp_ubar = 0;
 % Number of Load Increments i.e. Load is applied in "nsteps" increments 
loadingType = 'tension';
nsteps = 30;
forcevdisp = zeros(2,nsteps+1);
forcevdisp(1,1) = 0;
forcevdisp(2,1) = 0;

for step = 1 : nsteps 
        err3 = 1; %this is for intial tolerance.
        nit = 0;  %current iterations for this step
        Fint = zeros(total_unknown,1);
        fprintf(1,'\n Step %f \n',step);

        % -------------Original Montonic Loading-----------------
        % if step <= 5 
        %       % ubar = 0.012; %displacement for the steps <= 5  
        %       ubar = 1e-3;
        % elseif step > 5 
        %       ubar = 3e-3; %displacement for the steps > 5.
        % end
        ubar = 1e-3;
        temp_ubar = temp_ubar + ubar;
        disp([num2str(step)," ",num2str(temp_ubar)]);
    
        %Iterate until either the answer is converged or max iterations are 
        %reached.
        while ((err3>tol) && (nit<maxit))          % Newton Raphson loop
           
            nit = nit + 1; %incrementing the iterations
            
            % Computing Stiffness Matrix and Internal Force Vector
            [K,FAI,FE,D_st,kappa,NE_gp,stress_gp,interaction,stress_gp_sm,eq_stress,neq_stress,strain_gp_local,strain_gp_non_local,eq_strain] = peerlings_globalstifness(u_tot,strain_tot,D_st,material_p,De,damage_p,numelem,total_disp,...
                total_strain,node1,element1,element2,node2,elemType1,elemType2,kappa0,kappa,NE_gp,stress_gp);                     
                 
            Fint = [FAI;FE]; %internal forces vector
            R = Fint;
            
            % ------------------------------------------------
            %       Imposing Essential Boundary Conditions
            %-------------------------------------------------
    %         disp([num2str(toc),'   IMPOSING ESSENTIAL BOUNDARY CONDITION'])        
            
            Kres = K;
            bcwt = 1; % Used to keep the conditioning of the K matrix
            Kres(udofs,:) = 0;   % zero out the rows and columns of the K matrix, [As there is no displacement residual force is 0 ]
            Kres(:,udofs) = 0;   % zero out the corresponding columns             [ No displacement ]
            Kres(vdofs,:) = 0;   % making the rows and columns correspond to the nodes in lower boundary where displacement along Y-direction is not allowed
            Kres(:,vdofs) = 0;   
            Kres(udofs,udofs) = bcwt*speye(length(udofs));  % For inverse property
            Kres(vdofs,vdofs) = bcwt*speye(length(vdofs));  %   

            
            Kres(upper_disp_node,:) = 0; %There is no resdiual force for the vertical displacement of upper nodes, so corresponding columns are zero.
            Kres(:,upper_disp_node) = 0; % ??
            Kres(upper_disp_node,upper_disp_node) = bcwt*speye(length(upper_disp_node)); % For inverse property
            
            %In increment step,the displacement is applied to the upper_nodes so they get displaced by ubar. The equivalent extrenal force for this displacement is K(kk,zz)*ubar, which
            %should be equal to internal force, because we are sure the nodes are displaced by ubar and they experince internal force equivalent to external force for this displacement. 
            % So for the first iteration of increment change the residual.
            if nit == 1 
                for kk = 1:total_unknown               
                    for ll = 1:length(upper_disp_node)
                        zz = upper_disp_node(ll);
                        R(kk,:) = R(kk,:) - K(kk,zz)*ubar;
                    end
                end
            end
            
            R(udofs) = 0 ; % Imposing B.C. on Residual, As there is no displacement corresponding forces should also be zero.
            R(vdofs) = 0 ; % Imposing B.C. on Residual 
            
            if nit == 1
                R(upper_disp_node) = ubar;  % 
            else
                R(upper_disp_node) = 0;
            end 
            
            %------Solve for the correction---------%
            [du]= Kres\R; % Kres = 32805 x 32805 and R = 32805 x1 and du = 32805 x 1 [ Matrix right division ]
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
            STRAIN_LOCAL_XX(:,step) = strain_gp_local(:,1);  
            STRAIN_LOCAL_YY(:,step) = strain_gp_local(:,2);    
            STRAIN_LOCAL_XY(:,step)  = strain_gp_local(:,3);    
            STRAIN_NON_LOCAL_XX(:,step) = strain_gp_non_local(:,1); 
            STRAIN_NON_LOCAL_YY(:,step) = strain_gp_non_local(:,2);  
            STRAIN_NON_LOCAL_XY(:,step) = strain_gp_non_local(:,3); 
            EQ_STRAIN(:,step) = eq_strain;
            forcevdisp(1,step+1) = mean(u_tot(upper_disp_node,:));
            forcevdisp(2,step+1) = sum(Fint(lower_nodes,:));   
            forcevdisp(3,step+1) = sum(Fint(l_nodes,:)); 
            % 
            % save('Mode_I_steps1_80by80_Eta_4_R04_SmallLenScale_Beta9.mat','DAMAGE_DATA','NESTRAIN_DATA','GPT_DATA','forcevdisp','DISP_DATA','NESTRAIN_DATA_NODES','INTERNAL_FORCE',...
            %     'INTERACTION_DATA','SIGMA_XX','SIGMA_YY','SIGMA_XY','SIGMA_XX_smooth','SIGMA_YY_smooth','SIGMA_XY_smooth','EQ_STRESS','NEQ_STRESS');

             save(sprintf('Mode_I_steps_%d_%d_by_%d_Eta_%d_R04_SmallLenScale_Beta_%d.mat',nsteps,numx,numy,eta,beta),'DAMAGE_DATA','NESTRAIN_DATA','GPT_DATA','forcevdisp','DISP_DATA','NESTRAIN_DATA_NODES','INTERNAL_FORCE',...
                'INTERACTION_DATA','SIGMA_XX','SIGMA_YY','SIGMA_XY','SIGMA_XX_smooth','SIGMA_YY_smooth','SIGMA_XY_smooth','EQ_STRESS','NEQ_STRESS','node1','element1','STRAIN_LOCAL_XX', ...
                'STRAIN_LOCAL_YY','STRAIN_LOCAL_XY','STRAIN_NON_LOCAL_XX','STRAIN_NON_LOCAL_YY','STRAIN_NON_LOCAL_XY','EQ_STRAIN');
end
% end 
disp([num2str(toc),'  END OF NEWTON RAPHSON LOOP'])

% End of MaiN PrograM Code Function.
