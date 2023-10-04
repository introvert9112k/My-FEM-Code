
function [Gpnt] = computing_gauss_location(numelem,elemType,node2,element2)
gpnt = 0;
Gpnt = zeros(4*numelem,2);
Gpnt_Wt_detJ = zeros(4*numelem,1);
for iel = 1:numelem % Loop on elements
    order = 2;
    [W,Q] = quadrature(order,'GAUSS',2);  
    sctr2 = element2(iel,:); % Element connectivity 
    
  for kk = 1 : size(W,1) 
        gpnt = gpnt + 1;
        pt = Q(kk,:);   % Quadrature point 
        [N,dNdxi] = lagrange_basis(elemType,pt);
        J0 = node2(sctr2,:)'*dNdxi;                 % element Jacobian matrix       
        gpt_loc = N' * node2(sctr2,:);
        Gpnt(gpnt,:) = gpt_loc;  % GP in global coord, used
        Gpnt_Wt_detJ(gpnt,1) = W(kk);
        Gpnt_Wt_detJ(gpnt,2) = det(J0);
  end
end
end