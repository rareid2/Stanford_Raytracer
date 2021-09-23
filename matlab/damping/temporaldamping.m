function [gamma] = temporaldamping(f, kperp, kpar, w, m, wch, qh, mh, ...
                                   qs, Ns, ms, nus, B0, TOL, INTEGMETHOD)
% Compute the temporal damping rate gamma in a hot plasma for the
% resonances given in the vector m (m=0 is Landau damping, +-1 is the
% cyclotron resonance.
%
% f = cell array of distribution functions, each a function handle that
%     takes as input (vperp,vpar) and returns the value of the distribution
%     at that point.
% kperp:  perpendicular component of k vector
% kpar:   parallel component of k vector
% w:      frequency
% m:      resonances, a vector of resonances to sum in this integrand
% wch:    vector of hot (signed) gyrofrequencies, negative for electrons
% qh:     vector of hot charges, negative for electrons
% mh:     vector of hot masses
% qs:     vector of (cold) charges
% Ns:     vector of (cold) number densities in m^-3
% ms:     vector of (cold) masses
% nus:    vector of (cold) collision frequencies
% B0:     the magnetic field magnitude
% TOL:    integration tolerance
% INTEGMETHOD: integration method, either 'fast' (matlab builtin), or 
%              'accurate', a better quadrature method that almost always
%              converges, but is slower

DEL=1e-8;
d = DEL*abs(w);
if( d < 10*eps )
  d = 10*eps;
end;

dD0dw = ( hot_dispersion_real(kperp, kpar, w+d, qs, Ns, ms, nus, B0 ) - ...
          hot_dispersion_real(kperp, kpar, w-d, qs, Ns, ms, nus, B0 ) ) / (2*d);

gamma = 0;
for( ii=1:length(f) )
  Di = hot_dispersion_imag(f{ii}, kperp, kpar, w, m, wch(ii), ...
                           qh(ii), mh(ii), qs, Ns, ms, nus, B0, TOL, ... 
                           INTEGMETHOD );
  gamma = gamma + -Di/dD0dw;
end;
