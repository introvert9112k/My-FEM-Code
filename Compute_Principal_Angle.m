function [alpha] = Compute_Principal_Angle(stress_gpt)
    sigxx = stress_gpt(1,1);
    sigyy = stress_gpt(2,1);
    sigxy = stress_gpt(3,1);
    if sigxy==0 %already in principal orientation
        alpha = 0;
    else
        p = (2*sigxy)/(sigxx - sigyy); 
        alpha = 0.5*atand(p); %direct formula for calculating the principal angle.
    end
end
