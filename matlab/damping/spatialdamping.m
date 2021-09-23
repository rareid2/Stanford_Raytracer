function [ki] = spatialdamping(f, kperp, kpar, w, m, wch, qh, mh, ...
                               qs, Ns, ms, nus, B0, TOL, INTEGMETHOD )
% Compute the spatial damping rate ki in a hot plasma for the resonances
% given in the vector m (m=0 is Landau damping, +-1 is the cyclotron
% resonance.
%
% f = cell array of distribution functions, each a function handle that
%     takes as input (vperp,vpar) and returns the value of the distribution
%     at that point.
% kperp:  perpendicular component of k vector
% kpar:   parallel component of k vector
% w:      frequency
% m:      resonances, a vector of resonances to sum in this integrand
% qh:     vector of hot charges, negative for electrons
% mh:     vector of hot masses
% qs:     vector of charges
% Ns:     vector of number densities in m^-3
% ms:     vector of masses
% nus:    vector of collision frequencies
% B0:     the magnetic field magnitude
% TOL:    integration tolerance
% INTEGMETHOD: integration method, either 'fast' (matlab builtin), or 
%              'accurate', a better quadrature method that almost always
%              converges, but is slower

physconst;

% angle of wavenormal with respect to B0
theta = atan2(kperp,kpar);
% Refractive index
n = sqrt((clight^2/w^2)*(kperp^2+kpar^2));

[S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0);

A = S*sin(theta)^2 + P*cos(theta)^2;
B = R*L*sin(theta)^2 + P*S*(1+cos(theta)^2);

ki = 0;
for( ii=1:length(f) )
  Di = hot_dispersion_imag(f{ii}, kperp, kpar, w, m, wch(ii), ...
                           qh(ii), mh(ii), qs, Ns, ms, nus, B0, TOL, ...
                           INTEGMETHOD );
  ki = ki + -(w/clight)*(1/2)*(1/(4*n*(2*A*n^2-B))) * Di;
end;
