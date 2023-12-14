function U = element_disp(e,u,node1,element1)

% From the unknowns vector u, extract the parameters
% associated with the element "e"
% Then epsilon = B*U

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

end % END OF FUNCTION element_disp 