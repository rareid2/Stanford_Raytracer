function [k1,k2]=solve_dispersion_relation(n, w, qs, Ns, ms, nus, B0)
% Solve the dispersion relation and return k.  Inputs:
% 
% n = wavenormal direction
% w = frequency
% qs, Ns, ms, and nus are the charge, number density (m^-3), mass
% (kg), and collision frequency as column vectors, one per species.  B0 
% is the vector background magnetic field.
% 
% k1 and k2 are the positive and negative roots, respectively.

physconst;

k1 = zeros(size(w));
k2 = zeros(size(w));
for ii=1:length(w)
  % Find the angle the n vector makes with B0
  % cos^2(phi)
  cos2phi = ((n'*B0)*(n'*B0))/((n'*n)*(B0'*B0));
  sin2phi = 1.0 - cos2phi;
  
  % Magnitude of B0
  B0mag = norm(B0);
  
  % Find the stix parameters
  [S,D,P,R,L] = stix_parameters(w(ii), qs, Ns, ms, nus, B0mag);
  
  A = S*sin2phi+P*cos2phi;
  B = R*L*sin2phi+P*S*(1+cos2phi);
  nsquared1 = (B+sqrt(B^2-4*A*R*L*P))/(2*A);
  nsquared2 = (B-sqrt(B^2-4*A*R*L*P))/(2*A);
  
  n1 = sqrt(nsquared1);
  n2 = sqrt(nsquared2);
  
  k1(ii) = w(ii)*n1/c;
  k2(ii) = w(ii)*n2/c;
end