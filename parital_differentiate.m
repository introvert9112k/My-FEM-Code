clear
syms sigxx sigxy sigyy  
a = (sigxx + sigyy)/2;
b = (sigxx - sigyy)/2;
sig1 = a + sqrt((b^2) + (sigxy^2));
sig2 = a - sqrt((b^2) + (sigxy^2));
c11 = (sigxx^2 + sigxy^2);
c22 = (sigyy^2 + sigxy^2);
c12 = sigxy*(sigxx + sigyy);
Cbar = [c11 c12;
          c12 c22]/(sig2^2); 
d = diff(Cbar,sigxy)

% A1 = [sig1; sig2];
% A2 = [sigxx sigyy sigxy];
% d = diff(sig2,sigxy)

% dsig1_dsigxx = (sigxx/2 - sigyy/2)/(2*((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2)) + 1/2
% dsig1_dsigyy = 1/2 - (sigxx/2 - sigyy/2)/(2*((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2))
% dsig1_dsigxy = sigxy/((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2)

% dsig2_dsigxx = 1/2 - (sigxx/2 - sigyy/2)/(2*((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2))
% dsig2_dsigyy = (sigxx/2 - sigyy/2)/(2*((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2)) + 1/2
% dsig2_dsigxy = -sigxy/((sigxx/2 - sigyy/2)^2 + sigxy^2)^(1/2)
% sigxx = 0.0001;
% sigyy = 0.0002;
% sigxy = 0.0002;
% QQ = ((sigxx/2 - sigyy/2)^2 + sigxy^2);
% if QQ<0.0001
%     QQ = 0.0001;
% end
% aaa = (sigxx + sigyy)/(sigxx/2 + sigyy/2 - QQ^(1/2))^2 + (2*sigxy^2*(sigxx + sigyy))/(QQ^(1/2)*(sigxx/2 + sigyy/2 - QQ^(1/2))^3)

    