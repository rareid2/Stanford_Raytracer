function [sigma] = plasma_conductivity_tensor_symbolic(N,q,m,B0,nu,scale)
  physconst

  wp2 = N*q^2/m/EPS0;
  wc = q*B0/m;
  
  s = tf('s');
  if( nargin == 6 )
    s = s*scale;
  end
  sxx=s^2 + 2*nu*s + nu^2+wc(1)^2;
  sxy=wc(3)*s + wc(1)*wc(2)+nu*wc(3);
  sxz=-wc(2)*s - nu*wc(2)+wc(1)*wc(3);
  syx=-wc(3)*s + wc(1)*wc(2) - nu*wc(3);
  syy=s^2 + 2*nu*s + nu^2+wc(2)^2;
  syz=wc(1)*s + nu*wc(1)+wc(2)*wc(3);
  szx=wc(2)*s + nu*wc(2)+wc(1)*wc(3);
  szy=-wc(1)*s -nu*wc(1)+wc(2)*wc(3);
  szz=s^2 + 2*nu*s + nu^2+wc(3)^2;
  d = (s+nu)*(s^2+2*s*nu + nu^2+norm(wc)^2);

  sigma = EPS0*wp2/d * [sxx sxy sxz;
                        syx syy syz;
                        szx szy szz];
  