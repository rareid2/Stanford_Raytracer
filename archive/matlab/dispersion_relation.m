function [F]=dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
% Evaluate the dispersion relation function F(n,w) given the plasma
% parameters and the wavenormal n (in cartesian coordinates)
%
% n = refractive index vector
% w = frequency
% qs = vector of charges
% Ns = vector of number densities in m^-3
% ms = vector of masses
% nus = vector of collision frequencies
% B0 = the magnetic field (vector)

physconst

% Find the needed spherical components
%nmag2 = n'*n;
nmag2 = n(1,:).*n(1,:) + n(2,:).*n(2,:) + n(3,:).*n(3,:);
% cos^2(phi)
%cos2phi = ((n'*B0)*(n'*B0))/((n'*n)*(B0'*B0));
cos2phi = ( (n(1,:)*B0(1) + n(2,:)*B0(2) + n(3,:)*B0(3)).* ...
            (n(1,:)*B0(1) + n(2,:)*B0(2) + n(3,:)*B0(3)) ) ./ ...
          ( (n(1,:).*n(1,:) + n(2,:).*n(2,:) + n(3,:).*n(3,:))* ...
            (B0'*B0) );
sin2phi = 1.0 - cos2phi;

% Find the stix parameters
[S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, norm(B0));

% Old code
A = S*sin2phi + P*cos2phi;
B = R*L*sin2phi + P*S*(1+cos2phi);
F = A.*nmag2.^2 - B.*nmag2 + R*L*P;

