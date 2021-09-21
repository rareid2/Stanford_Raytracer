% phi,theta = wavenormal phi,theta
% B0 = vector background magnetic field
% Returns k
function [E,H, k] = plasma_modes_k(w, qs, Ns, ms, nus, B0, phi, theta)
  B0mag = norm(B0);
  
  physconst;

  [S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0mag);
  % Compute the unrotated permittivity tensor (aligned along z)
  perm0 = [S -j*D 0; j*D S 0; 0 0 P];
  % Rotation matrix
  [R] = rotation_matrix_z(B0);
  % Base permittivity 
  perm0 = R*perm0*R.';
  
  e = EPS0*perm0;
  u = MU0*eye(3);

  Khat = [0,                   -cos(phi),            sin(theta)*sin(phi);
          cos(phi),             0,                   -cos(theta)*sin(phi);
          -sin(theta)*sin(phi), cos(theta)*sin(phi), 0];
  
  % Form the matrix A
  Z = zeros(3,3);
  A = [Z, w*u; -w*e, Z];
  % Now B
  B = [Khat, Z; Z, Khat];
  
  % Solve for k and the modes
  [q,v] = eig(A,B);
  v = diag(v);

  [vsorted,sorti] = sort(abs(v));
  
  E = q(1:3,sorti(1:4));
  H = q(4:6,sorti(1:4));

  k = v(sorti(1:4));
  