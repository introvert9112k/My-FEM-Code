function plot_non_local_eq_strain_along_length 
%This function is for plotting the Values of Non local Equivalent Strain
%along the length i.e at each guass pont location along horizantal
%direction at specified load steps..
close all; clc;
% load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Tension.mat');
load('Mode_I_80by80_Eta_4_R04_SmallLenScale_Beta9.mat');
%load('Mode_I_steps_10_80_by_80_Eta_4_R04_SmallLenScale_Beta_9_Compression.mat');

L = 60; % Length of the plate
D = 60; % Width of the plate
numx = 80; % Number of elements in X direction
numy = 80; % Number of elements in Y direction
% Four corner points
nnx = numx+1; % Number of nodes in X-direction
nny = numy+1; % Number of nodes in Y-direction
pt1 = [0 0] ; pt2 = [L 0] ; pt3 = [L D] ; pt4 = [0 D] ;
elemType = 'Q4' ;
[node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy,elemType);

% ------------------For Q4 Elements-------------------%
numelem = size(element,1);

loadSteps = [30];

check_elem = ((numx/2)+1):numx; %we are intersted in only half section.
%This gives the element numbers in the second half.


figure
hold on 
%--------------Plot Non Local Equivalent Strain Variation----------------%
%  Showing the guass point numbering in the element.
%   -----------
%   3       1
%   4       2
%   -----------

%Plotting Non local equivalent strain variation along the length at Different load Steps
maxNonEquivalentStrain = 0;
finalMaxValue = 0;
for step = 1 : length(loadSteps)  
    it = 1;
    for del=1:1:numelem    
        if (ismember(del,check_elem))  %if element is in the second half        
                a = (del-1)*4;
                nonLocalEquivalentStrain(it,1) = NESTRAIN_DATA(a+4,loadSteps(step)); %getting the Damage value at that 4th guass point in element
                maxNonEquivalentStrain = max(maxNonEquivalentStrain,nonLocalEquivalentStrain(it,1));
                xcord(it,1) = GPT_DATA(a+4,1);   %x coordinate of the intersted guass point. 
                it=it+1;
        end 
    end
    finalMaxValue = max(finalMaxValue,maxNonEquivalentStrain);
    color = rand(1,3); %generates random rgb color
    plot((xcord/L),nonLocalEquivalentStrain,'Color' ,color,'LineStyle', '--','LineWidth',1);
end 

legends = cell(1,length(loadSteps));
for i = 1 : length(loadSteps) 
    legends{i} = num2str(loadSteps(i));
end 
axis([0.5 1 0 finalMaxValue + 0.05]) 
xlabel('x/L');
ylabel('Non Local Equivalent Strain');
title('Non Local Equivalent Strain along length at various load Steps');
legend(legends,'Location','best');
hold off
end 

