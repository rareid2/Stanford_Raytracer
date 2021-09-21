addpath('../../bin');
addpath('../');
TOL = 1e-3;
INTEGMETHOD = 'accurate';
DEL = 1e-8;


const;

% Background magnetic field
B0 = 2.4864e-7;

% Cold plasma density
Ne = 5e6;
% electron mass
me = ME;
% Temperature in eV 
kT = 5e3*Q; % 5 keV
% electron charge
qe = -Q;

% maxwellian distribution
Ne_h = .1e6;
fe = @(vperp,vpar) Ne_h*maxwellboltzmann(vperp,vpar,me,kT);
% Suprathermal distribution from Bell 2002
%fe = @(vperp,vpar) suprathermal(vperp,vpar);
% Hot cyclotron frequency
wce_h = (qe*B0/me);
% Hot charge
qe_h = qe;
% Hot mass
me_h = me;

Np = 5e6;
% proton mass
mp = MP;
% proton charge
qp = Q;

% vectors of species (cold plasma component)
qs = [qe qp];
Ns = [Ne Np];
nus = [0 0];
ms = [me mp];

% vector of hot plasma properties
qhs = [qe_h];
mhs = [me_h];
wchs = [wce_h];

% Resonances (0=landau)
m = [0];

wcs = ((qs*B0)./ms);
wps2 = (Ns.*qs.^2./ms./EPS0);
wps = sqrt(wps2);
    
wci=wcs(2);
wce=wcs(1);
wpi=wps(2);
wpe=wps(1);
wlhr = sqrt((power(wce,2) + power(wci,2) + power(wpe,2) + ...
             power(wpi,2) - sqrt(power(power(wce,2) - power(wci,2) + ...
                                       power(wpe,2),2) + ...
                                 2*(-power(wce,2) + power(wci,2) + ...
                                    power(wpe,2))*power(wpi,2) + ...
                                 power(wpi,4)))/2.);

ws = linspace(400,.7*abs(wcs(1)),100);

% cell array of distribution functions
fs = {fe};
  

plotnum=1;
for theta=2*pi/360*[0 10 20 30 40 50 60 70]

kis = zeros(size(ws));
for wi=1:length(ws)
  w = ws(wi);

  % Compute the stix parameters
  [S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0);
  
  % Solve the real dispersion relation to find kperp and kpar
  A = S*sin(theta)^2 + P*cos(theta)^2;
  B = R*L*sin(theta)^2 + P*S*(1+cos(theta)^2);
  C = R*L*P;
  n2 = roots([A,-B,C]);
  
  % Pick the whistler mode
  n = max(real(sqrt(n2)));
  
  fprintf('wi=%d of %d: n=%g\n', wi, length(ws), n);

  kmag = n*w/clight;
  
  if( kmag ~= 0  )
    kpar = kmag*cos(theta);
    kperp = kmag*sin(theta);

    % Compute the spatial damping rate
    ki =  spatialdamping( fs, kperp, kpar, w, m, wchs, qhs, mhs, ...
                          qs, Ns, ms, nus, B0, TOL, INTEGMETHOD );

    % COMPUTE THE GROUP VELOCITY
    d = DEL*abs(w);
    if( d < 10*eps )
      d = 10*eps;
    end;
    dD0dw = ( hot_dispersion_real(kperp, kpar, w+d, qs, ...
                                  Ns, ms, nus, B0 ) - ...
              hot_dispersion_real(kperp, kpar, w-d, qs, ...
                                  Ns, ms, nus, B0 ) ) / (2*d);
    d = DEL*abs(kpar);
    if( d < 10*eps )
      d = 10*eps;
    end;
    dD0dkpar = ( hot_dispersion_real(kperp, kpar+d, w, qs, ...
                                     Ns, ms, nus, B0 ) - ...
                 hot_dispersion_real(kperp, kpar-d, w, qs, ...
                                     Ns, ms, nus, B0 ) ) / (2*d);
    d = DEL*abs(kperp);
    if( d < 10*eps )
      d = 10*eps;
    end;
    dD0dkperp = ( hot_dispersion_real(kperp+d, kpar, w, qs, ...
                                      Ns, ms, nus, B0 ) - ...
                  hot_dispersion_real(kperp-d, kpar, w, qs, ...
                                      Ns, ms, nus, B0 ) ) / (2*d);
    
    vgpar = -dD0dkpar/dD0dw;
    vgperp = -dD0dkperp/dD0dw;
    
    krhatpar = kpar/sqrt(kpar^2+kperp^2);
    krhatperp = kperp/sqrt(kpar^2+kperp^2);
    vghatpar = vgpar/sqrt(vgpar^2+vgperp^2);
    vghatperp = vgperp/sqrt(vgpar^2+vgperp^2);

    kipar = ki*krhatpar;
    kiperp = ki*krhatperp;

    % Project ki along the group velocity vector
    kis(wi) = (kipar*vgpar+kiperp*vgperp)/sqrt(vgpar^2+vgperp^2);
  else
    kis(wi) = 0;
    disp('Re{n} = 0, not solving evanescent mode');
  end;

end;

wcs = ((qs*B0)./ms).*(w./(w+j*nus));
subplot(4,2,plotnum)
plot(abs(ws./wcs(1)), kis,'r')
axis tight
hold on
plotnum=plotnum+1;

end

set(gcf,'PaperSize',[10 7])
set(gcf,'PaperPosition',[0 0 10 7]);
print('-dpdf','blah.pdf')
