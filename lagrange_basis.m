function [Nv,dNdxi] = lagrange_basis(type,coord,dim)

    if ( nargin == 2 )
        dim=1;
    end

    switch type
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

end % % END OF FUNCTION lagrange_basis