function [f] = suprathermal(vperp,vpar)
  a = 4.9e5;
  b = 8.3e14;
  c = 5.4e23;

  % Just a crutch to avoid the singularity
  v0=1;

  % Convert to cm/s
  v = 100*sqrt(vperp.^2+vpar.^2+v0^2);
  f = (a./v.^4-b./v.^5+c./v.^6);

  % Convert to s^3/m^6 from s^3/cm^6 
  f = f*100^6;
  