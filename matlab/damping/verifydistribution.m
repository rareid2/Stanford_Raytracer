const
% Maxwellian
f = @(vperp,vpar) vperp.*maxwellboltzmann(vperp,vpar,ME,10e3*Q);
% Should be close to 1
blah=2*pi*dblquad(f,0,3e8,-3e8,3e8,1e-8)

