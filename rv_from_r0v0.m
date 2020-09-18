
function [R,V] = rv_from_r0v0(R0, V0, t)

global mu
%...Magnitudes of R0 and V0:
r0 = norm(R0);
v0 = norm(V0);
%...Initial radial velocity:
vr0 = dot(R0, V0)/r0;
%...Reciprocal of the semimajor axis (from the energy equation):
alpha = 2/r0 - v0^2/mu;
%...Compute the universal anomaly:
x = kepler_U(t, r0, vr0, alpha);
%...Compute the f and g functions:
[f, g] = f_and_g(x, t, r0, alpha);
%...Compute the final position vector:
R = f*R0 + g*V0;
%...Compute the magnitude of R:
r = norm(R);
%...Compute the derivatives of f and g:
[fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);
%...Compute the final velocity:
V = fdot*R0 + gdot*V0;
