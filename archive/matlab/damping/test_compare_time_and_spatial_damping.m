physconst;
addpath('../../bin');
addpath('../');

TOL = 1e-3;
INTEGMETHOD = 'accurate';

f = 400;
w = 2*pi*f;

% ray index
ri = 1;

out=readrayoutput('test_400Hz.txt');
t = out{ri}.time(1:10:end);

[pos,vprel,vgrel,n,B0,qs,ms,Ns,nus] = interpolaterayoutput(out{1},t);

% Ignore collisions for this code
nus = 0*nus;

% SET UP THE HOT ELECTRON DISTRIBUTION
% Hot electron density
Ne_h = .1e6;
% electron mass
me = ME;
% Temperature in eV 
kT = 5e3*Q; % 5 keV
% electron charge
qe = -Q;
% maxwellian distribution
fe = @(vperp,vpar) Ne_h*maxwellboltzmann(vperp,vpar,me,kT);
% Resonances (0=landau)
m = [0];
% cell array of hot distribution functions
fs = {fe};
  
% Test the spatial damping coefficient
magnitude_space = zeros(size(t));
magnitude_space(1) = 1;
% Test the time damping coefficient
magnitude_time = zeros(size(t));
magnitude_time(1) = 1;
for( ii=2:length(t) )
  % Hot cyclotron frequency
  wce_h = ((qe*norm(B0(ii,:)))./me);
  % Hot charge
  qe_h = qe;
  % Hot mass
  me_h = me;
  % vector of hot plasma properties
  wchs = [wce_h];
  qhs = [qe_h];
  mhs = [me_h];
  
  % Just trust that the n in the file is the one we want to use.
  % They might in fact be off a little due to the interpolation.
  k = n(ii,:)*w/clight;
  kmag = norm(k);
  % Signed component along B
  Bhat = B0(ii,:)/norm(B0(ii,:));
  kpar = k*Bhat';
  % Perpendicular vector
  kperp = k - kpar*Bhat;
  % Convert kperp into a magnitude
  kperp = norm(kperp);

  if( kmag ~= 0  )
    % Compute the spatial damping rate
    ki =  spatialdamping( fs, kperp, kpar, w, m, wchs, qhs, mhs, ...
                          qs(ii,:), Ns(ii,:), ms(ii,:), nus(ii,:), ...
                          norm(B0(ii,:)), TOL, INTEGMETHOD );
    gamma = temporaldamping( fs, kperp, kpar, w, m, wchs, qhs, mhs, ...
                             qs(ii,:), Ns(ii,:), ms(ii,:), nus(ii,:), ...
                             norm(B0(ii,:)), TOL, INTEGMETHOD );

    ki_along_vg = ki*(k*vgrel(ii,:)')/(norm(k)*norm(vgrel(ii,:)));

    dist = norm(pos(ii,:)-pos(ii-1,:));
    magnitude_space(ii) = magnitude_space(ii-1)*exp(-dist*ki_along_vg);
    magnitude_time(ii) = magnitude_time(ii-1)*exp(gamma*(t(ii)-t(ii-1)));
    fprintf('time=%g, space=%g\n', magnitude_time(ii), magnitude_space(ii));
  else
    disp('Re{n} = 0, not solving evanescent mode');
  end;
end;

