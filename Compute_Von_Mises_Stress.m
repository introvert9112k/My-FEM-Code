function [mises_gp] = Compute_Von_Mises_Stress(stress_gpt)
    sigxx = stress_gpt(1,1); 
    sigyy = stress_gpt(2,1);
    sigxy = stress_gpt(3,1);
    a = (sigxx + sigyy)/2;
    b = (sigxx - sigyy)/2;
    sig1 = a + sqrt((b^2) + (sigxy^2)); %obtained from the mohrs circle
    sig2 = a - sqrt((b^2) + (sigxy^2)); %obtained from the mohrs circle
    mises_gp = sqrt(0.5*((sig1-sig2)^2 + sig2^2 + sig1^2)); %von mises stress for principal plane stress case.
end 