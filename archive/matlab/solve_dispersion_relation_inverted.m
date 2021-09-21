function [w]=solve_dispersion_relation_inverted(n, k, qs, Ns, ms, nus, B0)
% Solve the dispersion relation and return w.  Inputs:
% 
% n = wavenormal direction
% k = spatial frequency
% qs, Ns, ms, and nus are the charge, number density (m^-3), mass
% (kg), and collision frequency as column vectors, one per species.  B0 
% is the vector background magnetic field.

physconst;

n = n/norm(n);
Khat = [  0   -n(3)   n(2);
         n(3)     0  -n(1);
        -n(2)  n(1)      0 ];

sz = 6+3*length(Ns);
w = zeros(length(k),sz);

matA = zeros(sz,sz);
matB = zeros(sz,sz);

% iterate over each species
offset = 7;
for jj=1:length(Ns)
  [sigma] = plasma_conductivity_tensor_symbolic(Ns(jj),qs(jj),ms(jj),B0,nus(jj));
  sigma=minreal(sigma);
  
  % Plasma conductivity
  sigma_ss = minreal(ss(sigma),[],false);
  A = sigma_ss.a;
  B = sigma_ss.b;
  C = sigma_ss.c;
  D = sigma_ss.d;
  
  matA(1:3,1:3) = matA(1:3,1:3)-D;
  matA(1:3,offset:(offset+2)) = -C;
  matA(offset:(offset+2),1:3) = B;
  matA(offset:(offset+2),offset:(offset+2)) = A;
  offset = offset+3;
end
matBinv = diag(1./[-i*EPS0*[1 1 1],i*MU0*[1 1 1],-i*ones(1,3*length(Ns))]);

for ii=1:length(k)
  matA(1:3,4:6) = i*k(ii)*Khat;
  matA(4:6,1:3) = i*k(ii)*Khat;
  
  matC = matBinv*matA;
  [v,d] = eig(matC);
  
  d=diag(d);
  [~,di]=sort((real(d)),'descend');
  d=d(di);
  v=v(:,di);
  w(ii,:) = d;
end
