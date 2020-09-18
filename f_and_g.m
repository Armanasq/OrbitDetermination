
function [f, g] = f_and_g(x, t, ro, a)

    %a - reciprocal of the semimajor axis (1/km)
    %ro - the radial position at time to (km)
    %t - the time elapsed since ro (s)
    %x - the universal anomaly after time t (km^0.5)
    %f - the Lagrange f coefficient (dimensionless)
    %g - the Lagrange g coefficient (s)

%%
    global mu
    z   = a*x^2;
    
%% Equ 3.69a
    f   = 1 - x^2/ro*stumpC(z);
    
%% Equ 3.69b
    g   = t - 1/sqrt(mu)*x^3*stumpS(z);
end
