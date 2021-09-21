function [ret] = integrand(f, vperp, kperp, kpar, w, m, wch, qh, mh, ...
                           R, L, P, S);
% f = distribution function, a function handle that takes as input 
%     (vperp,vpar) and returns the value of the distribution at that
%     point.
% vperp:  perpendicular velocity (may be a vector)
% kperp:  perpendicular component of k vector
% kpar:   parallel component of k vector
% w:      frequency
% m:      resonances, a vector of resonances to sum in this integrand
% wch:    hot particle (signed) gyrofrequency, negative for electrons
% qh:     hot particle (signed) charge, negative for electrons
% mh:     hot particle rest mass
% R:      stix R parameter (RH mode)
% L:      stix L parameter (LH mode)
% P:      stix P parameter (Parallel mode)
% S:      stix S parameter

physconst;

% angle of wavenormal with respect to B0
theta = atan2(kperp,kpar);
% Refractive index
n = sqrt((clight^2/w^2)*(kperp^2+kpar^2));

% cos(theta)
ct = cos(theta);
% sin(theta)
st = sin(theta);

vperp_vec = vperp;
ret = zeros(size(vperp_vec));
for( vperp_i=1:length(vperp_vec) )
  vperp = vperp_vec(vperp_i);

  % First sum the integrand over all m
  suminteg = 0;
  
  for mi=1:length(m)
    % J_{m}
    Jm = besselj(m(mi), kperp*vperp/wch);
    % J_{m-1}
    Jm_m1 = besselj(m(mi)-1, kperp*vperp/wch);
    % J_{m+1}
    Jm_p1 = besselj(m(mi)+1, kperp*vperp/wch);
    % This should of course be real, but matlab retruns a small imaginary
    % part sometimes even when it should not
    Jm = real(Jm);
    Jm_m1 = real(Jm_m1);
    Jm_p1 = real(Jm_p1);
    
    % Evaluate vpar only at the resonance
    vpar = (w-m(mi)*wch)/kpar;
    % evaluate G1
    G1 = fG1(f, vperp, vpar, kperp, kpar, w);
    % evaluate G2
    G2 = fG2(f, vperp, vpar, kperp, kpar, w, m(mi), wch);
    
    % evaluate the stuff under the sum

    % Chen's derivation.  Kennel has accidentally swapped (R-n^2) and (L-n^2)!
    suminteg = suminteg + ...
        (G1*((P-n^2*st^2)*(2*(L-n^2)*vperp*Jm_p1^2+2*vperp*(R-n^2)*Jm_m1^2+ ...
                           n^2*st^2*vperp*(Jm_p1-Jm_m1)^2) ...
             -n^2*ct*st*(2*vpar*Jm*(Jm_p1*(R-n^2)+Jm_m1*(L-n^2))+...
                         n^2*ct*st*vperp*(Jm_p1-Jm_m1)^2))...
         +G2*(4*vpar*Jm*((L-n^2)*(R-n^2)+n^2*st^2*(S-n^2))...
              -2*n^2*ct*st*((R-n^2)*vperp*Jm_m1+(L-n^2)*vperp*Jm_p1)));
  end;
  
  ret(vperp_i) = -2*pi^2*((qh^2/mh/EPS0)/(w*abs(kpar))) * suminteg*vperp;
end;
