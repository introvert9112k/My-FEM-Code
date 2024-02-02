function [K] = assemble_global_stiffness(total_unknown,sctr1,sctr3,Stif1,Stif2,Stif3,Stif4)
%-------------------- Assembling Global Stifness Matrix--------------
K = zeros(total_unknown,total_unknown);

    %--------------kaa---------------
    for i = 1 : size(sctr1,2)
        for j = 1 : size(sctr1,2) 
            ind1 = sctr1(1,i);
            ind2 = sctr1(1,j);
            K(ind1,ind2) = Stif1(i,j);
        end 
    end 
    
    %--------------kae---------------
    for i = 1 : size(sctr1,2)
        for j = 1 : size(sctr3,2) 
            ind1 = sctr1(1,i);
            ind2 = sctr3(1,j);
            K(ind1,ind2) = Stif2(i,j);
        end 
    end 

    %--------------kea---------------
    for i = 1 : size(sctr3,2)
        for j = 1 : size(sctr1,2) 
            ind1 = sctr3(1,i);
            ind2 = sctr1(1,j);
            K(ind1,ind2) = Stif3(i,j);
        end 
    end 
 
    %----------------kee---------------
    for i = 1 : size(sctr3,2)
        for j = 1 : size(sctr3,2) 
            ind1 = sctr3(1,i);
            ind2 = sctr3(1,j);
            K(ind1,ind2) = Stif4(i,j);
        end 
    end 

end 