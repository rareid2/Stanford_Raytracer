function [ret] = hot_dispersion_real(kperp, kpar, w, qs, Ns, ms, nus, B0)
% Compute the real component of the dispersion relation (same in 
% hot or cold)
%
% kperp:  perpendicular component of k vector
% kpar:   parallel component of k vector
% w:      frequency
% qs:     vector of charges
% Ns:     vector of number densities in m^-3
% ms:     vector of masses
% nus:    vector of collision frequencies
% B0:     the magnetic field magnitude

physconst;

[S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0);

% angle of wavenormal with respect to B0
theta = atan2(kperp,kpar);
% Refractive index
n = clight/w*sqrt(kperp^2+kpar^2);
A = S*sin(theta)^2 + P*cos(theta)^2;
B = R*L*sin(theta)^2 + P*S*(1+cos(theta)^2);
C = R*L*P;
ret = 4*(A*n^4 - B*n^2 + C);

% ret should equal the determinant of this.  It does.
% $$$ blah = [2*(R-n^2+1/2*n^2*sin(theta)^2), n^2*sin(theta)^2, n^2*cos(theta)* ...
% $$$         sin(theta);
% $$$         n^2*sin(theta)^2, 2*(L-n^2+1/2*n^2*sin(theta)^2), n^2*cos(theta)* ...
% $$$         sin(theta);
% $$$         n^2*cos(theta)*sin(theta), n^2*cos(theta)*sin(theta), P-n^2* ...
% $$$         sin(theta)^2];

