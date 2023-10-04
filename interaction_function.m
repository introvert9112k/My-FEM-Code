
function [g,dgdomega] = interaction_function(eta,omega,R)

    num = (1-R)*exp(-eta*omega) + R - exp(-eta);
    den = 1 - exp(-eta);
    g = num/den;
    num1 = (R-1)*eta*exp(-eta*omega);
    dgdomega = num1/den;

end % END OF FUNCTION interaction_function