% global stiffness
function [Stif,fai,fe,D_st,kappa,NE_gp,stress_gp,interaction,stress_gp_sm,eq_stress,neq_stress] = globalstiffness(u_tot,strain_tot,D_st,material_p, ...
    De,damage_p,numelem,total_disp,total_strain,node1,element1,element2,node2,elemType1,elemType2,kappa0,kappa,NE_gp,stress_gp) 

%-----------------Damage Parameters---------------------%
k = damage_p(1,1); 
alpha = damage_p(1,2);
beta = damage_p(1,3);
hcoup = damage_p(1,4);
R = damage_p(1,5);
eta = damage_p(1,6);
stress_gp_sm = zeros(4*numelem,3); %3 state of smoothed stress at gauss point
eq_stress = zeros(4*numelem,1); %equivalent stress at each gauss point
neq_stress = zeros(4*numelem,1); %non local equivalent stress at each gauss point.

%-----------------Material Parameters-------------------% 
nu = material_p(1,1); % Poisson Ratio 
E = material_p(1,2);
len_par = material_p(1,3); % Length Scale Parameter
th = material_p(1,4); % Thickness of the Specimen

del = [1; 1; 0];

H = [2 -1 0; % ####
    -1 2 0; 
     0 0 1.5];
 
tol = 1e-25; % Tolerance
tol1 = 1e-9;
tol2 = 1e-10;
tol3 = 1e-11;
total_unknown = total_disp + total_strain; % Total unknowns
gpnt = 0; % for storing the current guass point.
interaction = zeros(4*numelem,1); %interaction at each guass point location.
fai = zeros(total_disp,1); % Internal Force vector
fe = zeros(total_strain,1); % Internal Force vector corresponding to non-equivalent strains           ###
I = zeros(400*numelem,1); %  ###
J =  zeros(400*numelem,1); % ####
S = zeros(400*numelem,1); %  ###

index = 0; % ####

for iel = 1:numelem % Loop on elements  
    sctr = element2(iel,:); % location of the nodes of the current element in node1
    sctr1 = assembly(iel,element1);  %location of displacements of current element nodes in u_total [total displacements]
    sctr3 = assembly_nonlocal(sctr,total_disp);
    sctr2 = sctr3 - total_disp; %location of strains of current element nodes in strain_tot
    
    % Choose Gauss quadrature rules for elements
    order = 2;
    [W,Q] = quadrature(order,'GAUSS',2); %Q = quadrature points and W = quadrature weights

    u_local = u_tot(sctr1,:);  % Displacement corresponding to each node of the element  
    strain_nonlocal = strain_tot(sctr2,:); % Non-local Equivalent strain corresponding to each node of the element
    s1 = size(u_local,1); 
    s2 = size(strain_nonlocal,1);
    Stif1 = zeros(s1,s1); % 8*8
    Stif2 = zeros(s1,s2); % 8*12
    Stif3 = zeros(s2,s1); % 12*8
    Stif4 = zeros(s2,s2); % 12*12
    fai_gpt = zeros(s1,1); % Internal force at each guass point [2 components for each point ]
    fe_gpt = zeros(s2,1);  % 
    
    %iterating over each gauss point in element
    for kk = 1 : size(W,1) 
        
        gpnt = gpnt + 1;
        e_tilde1 = kappa(gpnt,1); % Non-equivalent strain at current gauss point of element   
        pt = Q(kk,:);   % Quadrature point       
        [B1,J01] = xfemBmatrix1(pt,elemType1,iel,node1,element1); % B matrix and Jacobian for displacement vector         
        [N2,B2] = xfemBmatrix2(pt,elemType2,iel,node2,element2); % Shape Function
        
        strain_gpt_local = B1*u_local; % Computing Strain at Gaussian Point
        strain_gpt_nonlocal = N2*strain_nonlocal; % Computing Non-local Strain vector at Quadrature point
        del_strain_gpt_nl = B2*strain_nonlocal; % Computing gradient of Non-local Strain vector at Quadrature point
        kappa0_gpt = kappa0(gpnt,1);      

        % Computing Nonlocal Equivalent Strain from Nonlocal Strain Vector
        % using modified von mises strain definition. 
        a1 = 1/(2*k);
        a2 = (k-1)/(1-2*nu);
        a3 = (2*k)/((1+nu)^2);
        a4 = 1/((1+nu)^2); 
        
        exx_l = strain_gpt_local(1,1);
        eyy_l = strain_gpt_local(2,1);
        exy_l = strain_gpt_local(3,1); 
        eyx_l = strain_gpt_local(3,1);
        I1_l = exx_l + eyy_l;
        J2_l = (2*(exx_l^2) + 2*(eyy_l^2) - 2*(exx_l*eyy_l)) + 3*(exy_l^2+eyx_l^2);   
        AA_l = (a2*I1_l)^2 + (a3)*J2_l;
        
        
        exx_nl = strain_gpt_nonlocal(1,1);
        eyy_nl = strain_gpt_nonlocal(2,1);
        exy_nl = strain_gpt_nonlocal(3,1); 
        eyx_nl = strain_gpt_nonlocal(3,1);
        I1_nl = exx_nl+eyy_nl;
        J2_nl = (2*(exx_nl^2) + 2*(eyy_nl^2) - 2*(exx_nl*eyy_nl)) + 3*(exy_nl^2+eyx_nl^2);   
        AA_nl = (a2*I1_nl)^2 + (a3)*J2_nl;
        
        if AA_l<=tol
            AA_l = tol;
        end
        
        if AA_nl<=tol
            AA_nl = tol;
        end
        
        eps_l = a1*(a2*I1_l + sqrt(AA_l)); % Computing Local Equivalent Strain from Strain Vector
        eps_nl = a1*(a2*I1_nl + sqrt(AA_nl)); % Computing Nonlocal Equivalent Strain from Strain Vector
        deps_de1 = a1*a2*(1 + a2*(AA_nl^(-0.5))*I1_nl)*del + a4*(AA_nl^(-0.5))*H*strain_gpt_nonlocal;
        deps_de_nl = deps_de1'; % Derivation of Nonlocal Equivalent strain wrt Nonlocal strain vector
        
        %----------------------Stress Smoothening to Compute Interaction Kernel----------------------------%   
        if eps_l == 0
            eps_l = tol;
        end
        
        if eps_nl == 0
            eps_nl = tol;
        end
        
        f = eps_nl - e_tilde1; % Damage Loading Function
        
        %if the non-local equivalent strain at current point > k,then
        %update k with current strain, because f <= 0.
        %kappa_gpt stores the maximum value of kappa achieved so far.
        if f>=0
            kappa_gpt = eps_nl;
        else
            kappa_gpt = e_tilde1;
        end
        
        % Exponential Cohesive Law i.e Damage Evolution Law         
        Omega = compute_damage(kappa_gpt,kappa0_gpt,alpha,beta);
        
        %------------------Constitutive Relations---------------------%
        %micromorphic stress vector.
        sm_stress_gpt = De*(1-Omega)*strain_gpt_nonlocal; % Smooth Stress Vector
%         sm_stress_gpt = De*strain_gpt_nonlocal; % Smooth Stress Vector
        
        if sm_stress_gpt== zeros(3,1)
            sm_stress_gpt(1,1) = tol1;
            sm_stress_gpt(2,1) = tol2;
            sm_stress_gpt(3,1) = tol3;
        end
        
        sigxx = sm_stress_gpt(1,1);
        sigyy = sm_stress_gpt(2,1);
        sigxy = sm_stress_gpt(3,1);                         
        
        %calculating micromorphic prinicpal stress components.
        [sigma1, sigma2, theta] = Compute_Principal_Stress(sm_stress_gpt);
        ptheta = (pi/180)*theta;
        Rot = [cos(ptheta) -sin(ptheta); sin(ptheta) cos(ptheta)];
                
        [g_int,dgdomega] = interaction_function(eta,Omega,R); % Computing Interaction Function 

        % Coefficient Omega
        if kappa_gpt < kappa0_gpt            
            DomegaDk = 0;  
        elseif e_tilde1 < eps_nl
            DomegaDk = (kappa0_gpt*(alpha*exp(beta*(kappa0_gpt - kappa_gpt)) - alpha + 1))/kappa_gpt^2 + (alpha*beta*kappa0_gpt*exp(beta*(kappa0_gpt - kappa_gpt)))/kappa_gpt;  
        else
            DomegaDk = 0;
        end    
        
        c_len = (len_par^2);
        
        QQ = ((sigxx/2 - sigyy/2)^2 + sigxy^2);
        
%         if QQ <= 0
%             QQ = 0.01;
%         end
              
        if (sigma1^2)>(sigma2^2)
        
        den = (sigma1^2);
        
        if den < 0.001
            den = 0.001;
        end
        
        c11 = (sigxx^2 + sigxy^2);
        c22 = (sigyy^2 + sigxy^2);
        c12 = (sigxy*(sigxx + sigyy));
        
        Cbar1 = [c11 c12; 
                  c12 c22]/den; 
              
        Cbar = [Cbar1 zeros(2,4);
                zeros(2,2) Cbar1  zeros(2,2);
                zeros(2,4) Cbar1];
               
        qq1 = [(2*sigxx) - (2*(sigxx^2 + sigxy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma1, sigxy - (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma1;
                sigxy - (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma1, -(2*(sigxy^2 + sigyy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma1];
        qq2 = [(2*(sigxx^2 + sigxy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma1,   sigxy + (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma1;
               sigxy + (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma1, (2*sigyy) + (2*(sigxy^2 + sigyy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma1];
        qq3 = [(2*sigxy) - (2*sigxy*(sigxx^2 + sigxy^2))/(QQ^(1/2)*sigma1), (sigxx + sigyy) - (2*sigxy^2*(sigxx + sigyy))/(QQ^(1/2)*sigma1);
               (sigxx + sigyy) - (2*sigxy^2*(sigxx + sigyy))/(QQ^(1/2)*sigma1),  (2*sigxy) - (2*sigxy*(sigxy^2 + sigyy^2))/(QQ^(1/2)*sigma1)]; 

        dCbar_dSigma(:,:,1) = [qq1 zeros(2,4);
                              zeros(2,2) qq1  zeros(2,2);
                              zeros(2,4) qq1]/den;
        dCbar_dSigma(:,:,2) = [qq2 zeros(2,4);
                               zeros(2,2) qq2  zeros(2,2);
                               zeros(2,4) qq2]/den;
        dCbar_dSigma(:,:,3) = [qq3 zeros(2,4);
                               zeros(2,2) qq3  zeros(2,2);
                               zeros(2,4) qq3]/den;                    
         
        dSigma_dStrain = (1-Omega)*De;  
        dSigma_dOmega = -(De*strain_gpt_nonlocal);
        Y2 = (dSigma_dStrain + dSigma_dOmega*DomegaDk*deps_de_nl)*N2;
        X1 = dCbar_dSigma;
        Y1 = del_strain_gpt_nl;
        E1_temp = cellfun(@(X1) X1*Y1, num2cell(X1,[1 2]),'UniformOutput',false);
        Eklp_1 = cat(3, E1_temp{:});
        Eklp = permute(Eklp_1,[1 3 2]);        
               
        dcdu(:,:) = Eklp*Y2;
        
        elseif (sigma2^2)>(sigma1^2)
        
        den = (sigma2^2);
        
        if den < 0.001
            den = 0.001;
        end
        
        c11 = (sigxx^2 + sigxy^2);
        c22 = (sigyy^2 + sigxy^2);
        c12 = (sigxy*(sigxx + sigyy));
        
        Cbar1 = [c11 c12;
                  c12 c22]/den; 

        Cbar =  [Cbar1 zeros(2,4);
                 zeros(2,2) Cbar1  zeros(2,2);
                 zeros(2,4) Cbar1];      
               
        qq1 = [(2*sigxx) + (2*(sigxx^2 + sigxy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma2, sigxy + (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma2;
                sigxy + (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma2,  (2*(sigxy^2 + sigyy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) - 1/2))/sigma2];
        qq2 =  [-(2*(sigxx^2 + sigxy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma2,   sigxy - (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma2;
               sigxy - (2*sigxy*(sigxx + sigyy)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma2, (2*sigyy) - (2*(sigxy^2 + sigyy^2)*((sigxx/2 - sigyy/2)/(2*QQ^(1/2)) + 1/2))/sigma2];
        qq3 =  [(2*sigxy) + (2*sigxy*(sigxx^2 + sigxy^2))/(QQ^(1/2)*sigma2), (sigxx + sigyy) + (2*sigxy^2*(sigxx + sigyy))/(QQ^(1/2)*sigma2);
               (sigxx + sigyy) + (2*sigxy^2*(sigxx + sigyy))/(QQ^(1/2)*sigma2),  (2*sigxy) + (2*sigxy*(sigxy^2 + sigyy^2))/(QQ^(1/2)*sigma2)]; 

        dCbar_dSigma(:,:,1) = [qq1 zeros(2,4);
                               zeros(2,2) qq1  zeros(2,2);
                               zeros(2,4) qq1]/den;
        dCbar_dSigma(:,:,2) = [qq2 zeros(2,4);
                               zeros(2,2) qq2  zeros(2,2);
                               zeros(2,4) qq2]/den;
        dCbar_dSigma(:,:,3) = [qq3 zeros(2,4);
                               zeros(2,2) qq3  zeros(2,2);
                               zeros(2,4) qq3]/den;                           
        
        dSigma_dStrain = (1-Omega)*De;  
        dSigma_dOmega = -(De*strain_gpt_nonlocal);
        Y2 = (dSigma_dStrain + dSigma_dOmega*DomegaDk*deps_de_nl)*N2;
        X1 = dCbar_dSigma;
        Y1 = del_strain_gpt_nl;
        E1_temp = cellfun(@(X1) X1*Y1, num2cell(X1,[1 2]),'UniformOutput',false);
        Eklp_1 = cat(3, E1_temp{:});
        Eklp = permute(Eklp_1,[1 3 2]);                      
        dcdu(:,:) = Eklp*Y2;
        
        end   
        
        C_grad = c_len*g_int*Cbar;
        
        % Constitutive Relations
        
        stress_gpt = (1-Omega)*De*strain_gpt_local + hcoup*(strain_gpt_local-strain_gpt_nonlocal);  % Stress corresponding to strain           
        
        conjugate_stress_gpt = hcoup*(strain_gpt_nonlocal - strain_gpt_local); % Conjugate Stress to Non-equivalent strain
        
        xit_gpt = hcoup*C_grad*del_strain_gpt_nl; % Conjugate Stress to Gradient of Non-equivalent strain
        
        
        % Computing Stiffness Matrices
        
        Coff1 = B1'*((1-Omega)*De + hcoup)*B1;
               
        Stif1 = Stif1 + Coff1*W(kk)*det(J01)*th;          
        
        Coff3 = -B1'*(hcoup + De*strain_gpt_local*DomegaDk*deps_de_nl)*N2; 
        
        Stif2 = Stif2 + Coff3*W(kk)*det(J01)*th; % Stiffness matrix Kaa (8 x 4)
        
        Coff5 = -(N2'*hcoup*B1);
               
        Stif3 = Stif3 + Coff5*W(kk)*det(J01)*th; % Stiffness matrix Kea (4 x 8);
        
        Coff6 = B2'*hcoup*C_grad*B2 + N2'*hcoup*N2 + B2'*hcoup*c_len*g_int*dcdu + B2'*hcoup*Cbar*del_strain_gpt_nl*dgdomega*DomegaDk*deps_de_nl*N2;
%         Coff6 = B2'*hcoup*C_grad*B2 + N2'*hcoup*N2 + B2'*hcoup*Cbar*del_strain_gpt_nl*dgdomega*DomegaDk*deps_de_nl*N2;

        Stif4 = Stif4 + (Coff6)*W(kk)*det(J01)*th; % Stiffness matrix Kea (4 x 4)  
        
        fai_gpt = fai_gpt - (B1'*stress_gpt)*W(kk)*det(J01)*th;
        
        fe_gpt = fe_gpt - (N2'*conjugate_stress_gpt + B2'*xit_gpt)*W(kk)*det(J01)*th;
        
        NE_gp(gpnt,1) = eps_nl;
        kappa(gpnt,1) = kappa_gpt;
        D_st(gpnt,1) = Omega;
        stress_gp(gpnt,:) = stress_gpt'; 
        interaction(gpnt,1) = g_int;       
        stress_gp(gpnt,:) = stress_gpt';
        stress_gp_sm(gpnt,:) = sm_stress_gpt';
        eq_stress(gpnt,:) = E*(1-Omega)*eps_l; 
        neq_stress(gpnt,:) = E*(1-Omega)*eps_nl;
    end                  % end of looping on GPs   
    
    stif_temp = [Stif1 Stif2; Stif3 Stif4];
    a = [sctr1 sctr3];
    b = [sctr1 sctr3];
    
    for ti = 1:20
        for tj = 1:20
            index = index + 1;
            I(index) =  a(ti);
            J(index) =  b(tj);
            S(index) = stif_temp(ti,tj);
        end
    end
            
%     Stif(a,b) = Stif(a,b) + stif_temp;
    
    fai(sctr1',1) = fai(sctr1',1) + fai_gpt;
    fe(sctr2',1) = fe(sctr2',1) + fe_gpt;
end                      % end of looping on elements

Stif = sparse(I,J,S,total_unknown,total_unknown);

end