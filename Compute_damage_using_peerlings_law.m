function [damage] = Compute_damage_using_peerlings_law(C,alpha,beta,neq_strain,D_i,kappa_gpt,kappa0_gpt)
if kappa_gpt <  kappa0_gpt 
     damage = 0;
else 
    damage = C*exp(alpha*D_i)*(neq_strain)^beta;
end 
end 