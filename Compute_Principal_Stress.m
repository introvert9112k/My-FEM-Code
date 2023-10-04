
function [sig1, sig2, alpha] = Compute_Principal_Stress(stress_gpt)
    sigxx = stress_gpt(1,1);
    sigyy = stress_gpt(2,1);
    sigxy = stress_gpt(3,1);
    a = (sigxx + sigyy)/2;
    b = (sigxx - sigyy)/2;
    sig1 = a + sqrt((b^2) + (sigxy^2));
    sig2 = a - sqrt((b^2) + (sigxy^2));
    den = (sigxx - sigyy);
    if den==0
        alpha = 90;
    else
        p = (2*sigxy)/den;
        alpha = 0.5*atand(p);
    end
end
