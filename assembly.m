function sctrB = assembly(e,element1)

sctr = element1(e,:); %gives the node numbers of the element.
nn   = length(sctr); % nn = 4

%sctrBfem is to store the location of the displaceemnt values for the each
%node of the element,becuase there are 2 displacements for each node.
for k = 1 : nn
    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
    sctrBfem(2*k)   = 2*sctr(k)   ;
end
    sctrB = sctrBfem;
end % End of  