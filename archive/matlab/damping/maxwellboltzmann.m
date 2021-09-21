function [f] = maxwellboltzmann(vperp,vpar,m,kT)
% Maxwell-boltzmann distribution in vperp,vpar space

f = (m/(2*pi*kT))^(3/2)*exp(-m*(vperp.^2+vpar.^2)/2/kT);

