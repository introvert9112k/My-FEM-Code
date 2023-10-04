%--------------------------Plot Output:Damage-----------------------------%
%-----------------Plot Output:Non-Equivalent Strain---------------------------%

function Compare_Load_Response_Small_Length

close all; clc;
L = 60; % Length of the plate
D = 60; % Width of the plate
step1 = [60 165 270];  % Small Length Scale Principle Stress Based Localizing GDM
step2 = [35 75 230]; % Small Length Scale Conventional Localizing GDM
loading = 'MODE_I';

%------------------Material Parameters------------------%
if (strcmp(loading,'MODE_I') )
    E  = 1000; % Elastic Moduli
    kappa0 = 0.002;
elseif (strcmp(loading,'COMPRESSION') )
    E  = 20000; % Elastic Moduli
    kappa0 = 0.0001;
end


%---------------------Load Displacement Plot---------------------%
figure
hold on
load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9_old','forcevdisp');
forcevdisp3 = forcevdisp;
plot((forcevdisp3(1,:)/L)*1e3,forcevdisp3(2,:)/(L*E*kappa0),'-r','LineWidth',1);
clear forcevdisp

load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9','forcevdisp');
forcevdisp1 = forcevdisp;
plot((forcevdisp1(1,:)/L)*1e3,forcevdisp1(2,:)/(L*E*kappa0),'-b','LineWidth',1);
clear forcevdisp
% plot(((forcevdisp1(1,step1(1,1)))/L)*1e3,(forcevdisp1(2,step1(1,1)))/(L*E*kappa0),'*k','LineWidth',1);
% plot(((forcevdisp1(1,step1(1,2)))/L)*1e3,(forcevdisp1(2,step1(1,2)))/(L*E*kappa0),'*k','LineWidth',1);
% plot(((forcevdisp1(1,step1(1,3)))/L)*1e3,(forcevdisp1(2,step1(1,3)))/(L*E*kappa0),'*k','LineWidth',1);
% plot(((forcevdisp3(1,step2(1,1)))/L)*1e3,(forcevdisp3(2,step2(1,1)))/(L*E*kappa0),'*k','LineWidth',1);
% plot(((forcevdisp3(1,step2(1,2)))/L)*1e3,(forcevdisp3(2,step2(1,2)))/(L*E*kappa0),'*k','LineWidth',1);
% plot(((forcevdisp3(1,step2(1,3)))/L)*1e3,(forcevdisp3(2,step2(1,3)))/(L*E*kappa0),'*k','LineWidth',1);
% load('Conventional_Mode_I_80by80_LargeLen_beta180_280818','forcevdisp');
% forcevdisp3 = forcevdisp;
% plot((forcevdisp3(1,:)/L)*1e3,forcevdisp3(2,:)/(L*E*kappa0),'--k','LineWidth',1);
% load('OUTPUT_ModeI_80by80_Beta_80_Eta_80_Large_Length_Scale_R05','forcevdisp');
% forcevdisp2 = forcevdisp;
% plot((forcevdisp2(1,:)/L)*1e3,forcevdisp2(2,:)/(L*E*kappa0),'-.k','LineWidth',1);
% clear forcevdisp
% legend('Proposed Model (l=0.04L)','Proposed Model (l=0.2L)','Conventional Model (l=0.04L)','Conventional Model (l=0.2L)','Location','northeast');
legend({'Conventional Model','Proposed Model'},'Location','northeast','FontSize',8);
% xlabel({'Displacement'},'FontSize',16);
% ylabel({'Force'},'FontSize',16);
axis([0 9 0 0.5]);
x0 = 100;
y0 = 100;
width=250;
height=175;
set(gcf,'units','points','position',[x0,y0,width,height])
end



function plot_mesh(X,connect,elem_type,se)

% function plot_mesh(X,connect,elem_type,linespec)
% plots a nodal mesh and an associated connectivity.  X is
% teh nodal coordinates, connect is the connectivity, and
% elem_type is either 'L2', 'L3', 'T3', 'T6', 'Q4', or 'Q9'
% depending on the element topology.

if ( nargin < 4 )
    se='w-';
end

holdState=ishold;
hold on

% fill X if needed
if (size(X,2) < 3)
    for c=size(X,2)+1:3
        X(:,c)=[zeros(size(X,1),1)];
    end
end

for e=1:size(connect,1)

    if ( strcmp(elem_type,'Q9') )       % 9-node quad element
        ord=[1,5,2,6,3,7,4,8,1];
    elseif ( strcmp(elem_type,'Q8') )  % 8-node quad element
        ord=[1,5,2,6,3,7,4,8,1];
    elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
        ord=[1,2,3,1];
    elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
        ord=[1,4,2,5,3,6,1];
    elseif ( strcmp(elem_type,'Q4') )  % 4-node quadrilateral element
        ord=[1,2,3,4,1];
    elseif ( strcmp(elem_type,'L2') )  % 2-node line element
        ord=[1,2];
    elseif ( strcmp(elem_type,'L3') )  % 3-node line element
        ord=[1,3,2];
    elseif ( strcmp(elem_type,'H4') )  % 4-node tet element
        ord=[1,2,4,1,3,4,2,3];
    elseif ( strcmp(elem_type,'B8') )  % 8-node brick element
        ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];
    end

    for n=1:size(ord,2)
        xpt(n)=X(connect(e,ord(n)),1);
        ypt(n)=X(connect(e,ord(n)),2);
        zpt(n)=X(connect(e,ord(n)),3);
    end
    plot3(xpt,ypt,zpt,se)
end

rotate3d on
axis equal

if ( ~holdState )
    hold off
end
end % END OF FUNCTION plotmesh



function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)

% function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)
%
% creates a connectivity list

if ( nargin < 5 )
    disp(['Not enough parameters specified for make_elem function'])
end

inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v,size(node_pattern,2));

for row=1:num_v
    for col=1:num_u
        element(e,:)=node_pattern+inc;
        inc=inc+inc_u;
        e=e+1;
    end
    inc=row*inc_v;
end
end

function sctrB = assembly_nonlocal(sctr,total_disp_local)

nn   = length(sctr); % Number of nodes
sctr_n = zeros(1,nn);
for i = 1:nn
a = sctr(i);    
b = a + total_disp_local;
sctr_n(i) = b + a*3;
end

sctr_nonlocal = zeros(1,16);

cnt = 0 ;

for k = 1 : nn
    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k) - 3;

    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k) - 2;

    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k) - 1;
    
    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k);
end
sctrB = sctr_nonlocal;
end

function sctrB = assembly(sctr)
sctrBfem = zeros(1,8);
nn   = length(sctr);
for k = 1 : nn
    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
    sctrBfem(2*k)   = 2*sctr(k)   ;
end
    sctrB = sctrBfem;
end

function U = element_disp(e,u)

% From the unknowns vector u, extract the parameters
% associated with the element "e"
% Then epsilon = B*U

global node1 element1

sctr = element1(e,:);
nn   = length(sctr);

% stdU contains true nodal displacement
idx = 0 ;
stdU   = zeros(2*nn,1);
for in = 1 : nn
    idx = idx + 1;
    nodeI = sctr(in) ;
    stdU(2*idx-1) = u(2*nodeI-1);
    stdU(2*idx)   = u(2*nodeI  );
end

U = stdU;
end

function [Nv,dNdxi]=lagrange_basis(type,coord,dim)

% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
%   H4, H10, B8 and B27
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.

if ( nargin == 2 )
    dim=1;
end

switch type
    case 'L2'
        %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
        %
        %    1---------2
        %
        if size(coord,2) < 1
            disp('Error coordinate needed for the L2 element')
        else
            xi=coord(1);
            N=([1-xi,1+xi]/2)';
            dNdxi=[-1;1]/2;
        end

    case 'L3'
        %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
        %
        %    1---------2----------3
        %
        if size(coord,2) < 1
            disp('Error two coordinates needed for the L3 element')
        else
            xi=coord(1);
            N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2];
            dNdxi=[xi-.5;xi+.5;-2*xi];
        end

    case 'T3'
        %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
        %
        %               3
        %             /  \
        %            /    \
        %           /      \
        %          /        \
        %         /          \
        %        /            \
        %       /              \
        %      /                \
        %     /                  \
        %    1--------------------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the T3 element')
        else
            xi=coord(1); eta=coord(2);
            N=[1-xi-eta;xi;eta];
            dNdxi=[-1,-1;1,0;0,1];
        end

    case 'T3fs'
        if size(coord,2) < 2
            disp('Error two coordinates needed for the T3fs element')
        else
            xi=coord(1); eta=coord(2);
            N=[1-xi-eta;xi;eta];
            dNdxi=[-1,-1;1,0;0,1];
        end

    case 'T4'
        %%%%%%%%%% T4 FOUR NODE TRIANGULAR CUBIC BUBBLE ELEMENT %%%%%%%%%%%%
        %
        %               3
        %             /  \
        %            /    \
        %           /      \
        %          /        \
        %         /          \
        %        /      4     \
        %       /              \
        %      /                \
        %     /                  \
        %    1--------------------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the T4 element')
        else
            xi=coord(1); eta=coord(2);
            N=[1-xi-eta-3*xi*eta;xi*(1-3*eta);eta*(1-3*xi);9*xi*eta];
            dNdxi=[-1-3*eta,-1-3*xi;
                1-3*eta, -3*xi;
                -3*eta,   1-3*xi;
                9*eta,   9*xi ];
        end

    case 'T6'
        %%%%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
        %
        %               3
        %             /  \
        %            /    \
        %           /      \
        %          /        \
        %         6          5
        %        /            \
        %       /              \
        %      /                \
        %     /                  \
        %    1---------4----------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the T6 element')
        else
            xi=coord(1); eta=coord(2);
            N=[1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                xi*(2*xi-1);
                eta*(2*eta-1);
                4*xi*(1-xi-eta);
                4*xi*eta;
                4*eta*(1-xi-eta)];

            dNdxi=[4*(xi+eta)-3   4*(xi+eta)-3;
                4*xi-1              0;
                0        4*eta-1;
                4*(1-eta-2*xi)          -4*xi;
                4*eta           4*xi;
                -4*eta  4*(1-xi-2*eta)];
        end


    case 'Q4'
        %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4--------------------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1--------------------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the Q4 element')
        else
            xi=coord(1); eta=coord(2);
            N=1/4*[ (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];
            dNdxi=1/4*[-(1-eta), -(1-xi);
                1-eta,    -(1+xi);
                1+eta,      1+xi;
                -(1+eta),   1-xi];
        end
      case 'Q8'
        %%%%%%%%%%%%%%% Q8 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4---------7----------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    8                    6
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1---------5----------2
        %
        if size(coord,2) < 2
        disp('Error two coordinates needed for the Q8 element')
        else
        xi=coord(1); eta=coord(2);
        N=1/4*[-1*(1-xi)*(1-eta)*(1+xi+eta);
                -1*(1+xi)*(1-eta)*(1-xi+eta);
                -1*(1+xi)*(1+eta)*(1-xi-eta);
                -1*(1-xi)*(1+eta)*(1+xi-eta);
                2*(1-xi^2)*(1-eta);
                2*(1+xi)*(1-eta^2);
                2*(1-xi^2)*(1+eta);
                2*(1-xi)*(1-eta^2)];
        dNdxi=1/4*[(1-eta)*(2*xi+eta),  (1-xi)*(2*eta+xi);
                    (1-eta)*(2*xi-eta), (1+xi)*(2*eta-xi);
                    (1+eta)*(2*xi+eta), (1+xi)*(2*eta+xi);
                    (1+eta)*(2*xi-eta), (1-xi)*(2*eta-xi);
                    -4*xi*(1-eta),      -2*(1-xi^2);
                    2*(1-eta^2),        -4*eta*(1+xi);
                    -4*xi*(1+eta),       2*(1-xi^2);
                    -2*(1-eta^2),       -4*eta*(1-xi)];
        end
    case 'Q9'
        %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4---------7----------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    8          9         6
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1----------5---------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the Q9 element')
        else
            xi=coord(1); eta=coord(2);
            N=1/4*[xi*eta*(xi-1)*(eta-1);
                xi*eta*(xi+1)*(eta-1);
                xi*eta*(xi+1)*(eta+1);
                xi*eta*(xi-1)*(eta+1);
                -2*eta*(xi+1)*(xi-1)*(eta-1);
                -2*xi*(xi+1)*(eta+1)*(eta-1);
                -2*eta*(xi+1)*(xi-1)*(eta+1);
                -2*xi*(xi-1)*(eta+1)*(eta-1);
                4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
            dNdxi=1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
                eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
                eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
                eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
                -4*xi*eta*(eta-1),   -2*(xi+1)*(xi-1)*(2*eta-1);
                -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
                -4*xi*eta*(eta+1),   -2*(xi+1)*(xi-1)*(2*eta+1);
                -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
                8*xi*(eta^2-1),      8*eta*(xi^2-1)];
        end

    case 'H4'
        %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
        %
        %             4
        %           / | \
        %          /  |  \
        %         /   |   \
        %        /    |    \
        %       /     |     \
        %      1 -----|------3
        %         -   2  -
        if size(coord,2) < 3
            disp('Error three coordinates needed for the H4 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            N=[1-xi-eta-zeta;
                xi;
                eta;
                zeta];
            dNdxi=[-1  -1  -1;
                1   0   0;
                0   1   0;
                0   0   1];
        end

    case 'H10'
        %%%%%%%%%%%%%%%% H10 TEN NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
        disp(['Element ',type,' not yet supported'])
        if size(coord,2) < 3
            disp('Error three coordinates needed for the H10 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            N=zeros(10,1);
            dNdxi=zeros(10,3);
        end

    case 'B8'
        %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
        %
        %                  8
        %               /    \
        %            /          \
        %         /                \
        %      5                     \
        %      |\                     7
        %      |   \                / |
        %      |     \     4    /     |
        %      |        \    /        |
        %      |           6          |
        %      1           |          |
        %       \          |          3
        %          \       |        /
        %            \     |     /
        %               \  |  /
        %                  2
        %
        if size(coord,2) < 3
            disp('Error three coordinates needed for the B8 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            I1=1/2-coord/2;
            I2=1/2+coord/2;
            N=[   I1(1)*I1(2)*I1(3);
                I2(1)*I1(2)*I1(3);
                I2(1)*I2(2)*I1(3);
                I1(1)*I2(2)*I1(3);
                I1(1)*I1(2)*I2(3);
                I2(1)*I1(2)*I2(3);
                I2(1)*I2(2)*I2(3);
                I1(1)*I2(2)*I2(3)   ];
            dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;
        end

    case 'B27'
        %%%%%%%%%%%%%% B27 TWENTY SEVEN NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%
        disp(['Element ',type,' not yet supported'])
        if size(coord,2) < 3
            disp('Error three coordinates needed for the B27 element')
        else
            xi=coord(1); eta=coord(2); zeta=coord(3);
            N=zeros(27,1);
            dNdxi=zeros(27,3);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        disp(['Element ',type,' not yet supported'])
        N=[]; dNdxi=[];
end

I=eye(dim);
Nv=[];
for i=1:size(N,1)
    Nv=[Nv;I*N(i)];
end

if ( dim == 1 )
    B=dNdxi;
elseif ( dim == 2 )
    B=zeros(dim*size(N,1),3);

    B(1:dim:dim*size(N,1)-1,1) = dNdxi(:,1);
    B(2:dim:dim*size(N,1),2)   = dNdxi(:,2);

    B(1:dim:dim*size(N,1)-1,3) = dNdxi(:,2);
    B(2:dim:dim*size(N,1),3)   = dNdxi(:,1);
elseif ( dim == 3 )
    B=zeros(dim*size(N,1),6);

    disp('Error: need to add 3D N and dNdxi')

    B(1:dim:dim*size(N,1)-2,1) = dNdxi(:,1);
    B(2:dim:dim*size(N,1)-1,2) = dNdxi(:,2);
    B(3:dim:dim*size(N,1),3)   = dNdxi(:,3);

    B(2:dim:dim*size(N,1)-1,4) = dNdxi(:,3);
    B(3:dim:dim*size(N,1),4)   = dNdxi(:,2);

    B(3:dim:dim*size(N,1),5)   = dNdxi(:,1);
    B(1:dim:dim*size(N,1)-2,5) = dNdxi(:,3);

    B(1:dim:dim*size(N,1)-2,6) = dNdxi(:,2);
    B(2:dim:dim*size(N,1)-1,6) = dNdxi(:,1);
end
end % end of function


function plot_field(X,connect,elem_type,field)

% function plot_field(X,connect,elem_type,field)

if ( nargin == 4 )
    nodesoff=0;
end

if ( size(field) == size(connect) )
    elementalField=1;
else
    elementalField=0;
end

% fill X if needed
if (size(X,2) < 3)
    for c=size(X,2)+1:3
        X(:,c)=[zeros(size(X,1),1)];
    end
end

holdState=ishold;
hold on

% plot elements
if ( strcmp(elem_type,'Q9') )      % Q9 element
    ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'T3') )  % T3 element
    ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T4') )  % T4 element
    ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T6') )  % T6 element
    ord=[1,4,2,5,3,6,1];
elseif ( strcmp(elem_type,'Q4') )  % Q4 element
    ord=[1,2,3,4,1];
elseif ( strcmp(elem_type,'Q8') )  % T3 element
    ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'L2') )  % L2 element
    ord=[1,2];
elseif ( strcmp(elem_type,'L3') )  % L3 element
    ord=[1,3,2];
end

for e=1:size(connect,1)

    xpt=X(connect(e,ord),1);
    ypt=X(connect(e,ord),2);
    zpt=X(connect(e,ord),3);

    if ( elementalField )
        fpt=field(e,ord);
    else
        fpt=field(connect(e,ord));
    end

    fill3(xpt,ypt,zpt,fpt)
end

shading interp
axis equal

if ( ~holdState )
%     hold off
end
end

function conn=tricheck(node,conn,verbose)

% FUNCTION
%
%   conn=tricheck(node,conn,verbose)
%
% This function check wether a triangle has a negative Jacobian, and if
% so reorders it so that the the Jacobian is positive.

if ( nargin==2 )
    verbose=0;
end

if ( size(node,2)==3 )
    node=node(:,1:2);
end

count=0;

for e=1:size(conn,1)

    sctr=conn(e,:);
    [N,dNdxi]=lagrange_basis('T3',[1/3 1/3]);
    detJ=det(node(sctr,:)'*dNdxi);

    if ( detJ < 0 )
        %disp(['NEGATIVE JACOBIAN IN ELEMENT ',num2str(e)])
        conn(e,:)=fliplr(sctr);
        count=count+1;
    elseif ( detJ == 0 )
        disp(['ZERO JACOBIAN IN ELEMENT ',num2str(e),' CANNOT FIX'])
    end
end

if ( verbose )
    disp(['TRICHECK FOUND ',num2str(count),' NEGATIVE JACOBIANS, ALL FIXED'])
end
end

function [node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy, elemType)

switch elemType

    case 'Q4'           % here we generate the mesh of Q4 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);

    case 'Q8'           % here we generate a mesh of Q9 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
        [element,node]=q4totq8(element,node,numx,numy);

    otherwise
        error('For now, only PARTICLE, Q4, Q9 and T3 are supported by the mesh generator');
end
end
function X=square_node_array(pt1,pt2,pt3,pt4,numnod_u,numnod_v,uratio,vratio)

% square_node_array
%
% Generates a quadratleral array of nodes between the counterclockwise
% ordering of nodes pt1 - pt4.  There are numnod_u nodes in the u direction
% (pt1 - pt2) and numnode_v nodes in the v direction (pt2 - pt3).  The
% parameter uratio and vratio determint the nodal spacing along the u and v
% lines.  If no values of uratio and/or vratio are given values of unity are
% assumed which resulets in uniformed node spacng along the u an v directions.
% uratio and v ratio  are the ratio of the first node spacing to the last
% node spacing along th u of v direction respectivly (the first spacing
% occurs near pt 1 and teh last near pt3.
%
% X=square_node_array(pt1,pt2,pt3,pt4,numnod_u,numnod_v,uratio,vratio)
%

if ( nargin < 6 )
    disp(['Not enough parameters specified for quare_node_array function'])
elseif ( nargin == 6 )
    uratio=1;
    vratio=1;
elseif ( nargin == 7 )
    vratio=1;
end

% get node spacing along u direction
if ( uratio == 1 )
    xi_pts=linspace(-1,1,numnod_u);
elseif ( uratio > 0 )
    ru=uratio^(1/(numnod_u-2));
    xi_pts(1)=0;
    d=1;
    for i=2:numnod_u
        xi_pts(i)=xi_pts(i-1)+d;
        d=d/ru;
    end
    xi_pts=2*xi_pts/xi_pts(numnod_u)-1;
else
    disp('uratio must be greator than 0');
    xi_pts=linspace(-1,1,numnod_u);
end

% get node spacing along v direction
if ( vratio == 1 )
    eta_pts=linspace(-1,1,numnod_v);
elseif ( vratio > 0 )
    rv=vratio^(1/(numnod_v-2));
    eta_pts(1)=0;
    d=1;
    for i=2:numnod_v
        eta_pts(i)=eta_pts(i-1)+d;
        d=d/rv;
    end
    eta_pts=2*eta_pts/eta_pts(numnod_v)-1;
else
    disp('vratio must be greator than 0');
    eta_pts=linspace(-1,1,numnod_v);
end

x_pts=[pt1(1),pt2(1),pt3(1),pt4(1)];
y_pts=[pt1(2),pt2(2),pt3(2),pt4(2)];

for r=1:numnod_v
    eta=eta_pts(r);
    for c=1:numnod_u
        xi=xi_pts(c);
        % get interpolation basis at xi, eta
        N=lagrange_basis('Q4',[xi,eta]);
        N=N(:,1);
        X((r-1)*numnod_u+c,:)=[x_pts*N,y_pts*N];
    end
end
end