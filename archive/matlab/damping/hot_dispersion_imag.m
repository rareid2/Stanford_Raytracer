function [ret] = hot_dispersion_imag(f, kperp, kpar, w, m, wch, qh, mh, ...
                                     qs, Ns, ms, nus, B0, TOL, INTEGMETHOD)
% Compute the imaginary component of the hot dispersion relation for the 
% resonances given in the vector m (m=0 is Landau damping, +-1 is the 
% cyclotron resonance.
%
% f = distribution function, a function handle that takes as input 
%     (vperp,vpar) and returns the value of the distribution at that
%     point.
% kperp:  perpendicular component of k vector
% kpar:   parallel component of k vector
% w:      frequency
% m:      resonances, a vector of resonances to sum in this integrand
%         0=landau resonance, -1,+1=cyclotron resonances
% wch:    hot (signed) gyrofrequency, negative for electrons
% qh:     hot species charge
% mh:     hot species mass
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

[S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0);

% adaptive integration of the integrand over vperp, [0, inf]
integrand_vperp = @(vperp) integrand( f, vperp, kperp, kpar, w, m, ...
                                      wch, qh, mh, ...
                                      R, L, P, S );

% Normalize by approximately c
k=clight;
integrand_vperpnorm = @(vperp) k*integrand_vperp(vperp*k);

% Remap into a finite range [0+,1]
integrand_t = @(t) ...
    ((1+eps)./(t.^2+eps)).*integrand_vperpnorm((1-t+eps)./(t+eps));

if( strcmp( INTEGMETHOD, 'accurate' ) )
  % This method is pretty reliable.
  [integrated_integrand,errbnd] = quadva(integrand_t,[0,1],TOL,eps);
elseif( strcmp( INTEGMETHOD, 'fast' ) )
  % Integrate
  % matlab builtin struggles with this function
  integrated_integrand = quadl(integrand_t,0,1,TOL);
% $$$   disp('test')
% $$$   integrated_integrand = quad(integrand_t,0,1,TOL)
% $$$   integrated_integrand = quadl(integrand_t,0,1,TOL)
% $$$   integrated_integrand = quadgk(integrand_t,0,1)
% $$$   integrated_integrand = quadva(integrand_t,[0,1],TOL,eps)
% $$$   t=linspace(0,1,10000);
% $$$   plot(t,integrand_t(t))
% $$$   pause
end;

%ret = -2*pi^2*((qh^2/mh/EPS0)/(w*abs(kpar))) * integrated_integrand;
ret = integrated_integrand;
