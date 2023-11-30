function plot_damage_variation_gauss_point_along_length_1
%--------------------------Plot Output:Damage-----------------------------%

close all; clc;
load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9.mat');
L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 80; % Number of elements in X direction
numy = 80; % Number of elements in Y direction
% % Four corner points
nnx = numx+1; % Number of nodes in X-direction
nny = numy+1; % Number of nodes in Y-direction
pt1 = [0 0] ; pt2 = [L 0] ; pt3 = [L D] ; pt4 = [0 D] ;
elemType = 'Q4' ;
[node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType);

%------------------Material Parameters------------------%
loading = 'MODE_I';
% loading = 'COMPRESSION';
if (strcmp(loading,'MODE_I') )
    E  = 1000; % Elastic Moduli
    kappa0 = 0.002;
elseif (strcmp(loading,'COMPRESSION') )
    E  = 20000; % Elastic Moduli
    kappa0 = 0.0001;
end

% % ------------------For Q4 Elements-------------------%

numelem = size(element,1);

loadSteps = [5,10,15,20,25,30];

check_elem = ((numx/2)+1):numx; %we are intersted in only half section.
%This gives the element numbers in the second half.

%---------------------Load Displacement Plot---------------------%

% subplot(1,2,1)

plot((forcevdisp(1,:)/L)*1e3,forcevdisp(2,:)/(L*E*kappa0),'--k','LineWidth',1);
axis([0 9 0 0.6]);
hold on
plot(((forcevdisp(1,loadSteps(1)))/L)*1e3,(forcevdisp(2,loadSteps(1)))/(L*E*kappa0),'*k','LineWidth',1);
hold on
plot(((forcevdisp(1,loadSteps(5)))/L)*1e3,(forcevdisp(2,loadSteps(5)))/(L*E*kappa0),'*k','LineWidth',1);
hold on
plot(((forcevdisp(1,loadSteps(6)))/L)*1e3,(forcevdisp(2,loadSteps(6)))/(L*E*kappa0),'*k','LineWidth',1);
hold off


figure
hold on 
%--------------Plot Damage Variation----------------%
%  Showing the guass point numbering in the element.
%   -----------
%   3       1
%   4       2
%   -----------

%Plotting damage variation along the length at Different load Steps
for step = 1 : length(loadSteps)  
    it = 1;
    for del=1:1:numelem    
        if (ismember(del,check_elem))  %if element is in the second half        
                a = (del-1)*4;
                omega(it,1) = DAMAGE_DATA(a+4,loadSteps(step)); %getting the Damage value at that 4th guass point in element
                xcord(it,1) = GPT_DATA(a+4,1);   %x coordinate of the intersted guass point. 
                it=it+1;
        end 
    end
    color = rand(1,3); %generates random rgb color
    plot((xcord/L),omega,'Color' ,color,'LineStyle', '--','LineWidth',1);
end 

legends = cell(1,length(loadSteps));
for i = 1 : length(loadSteps) 
    legends{i} = num2str(loadSteps(i));
end 
axis([0.5 1 0 1]) 
xlabel('x/L');
ylabel('Damage');
% title('Damage along length at various load Steps');
legend(legends,'Location','best');
hold off
end 

