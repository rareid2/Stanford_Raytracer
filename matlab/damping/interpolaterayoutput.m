function [pos,vprel,vgrel,n,B0,qs,ms,Ns,nus] = interpolaterayoutput(raystruct,t)
% Interpolate one output from readrayoutput over the given t (time) vector
% Interpolates everything:  position, etc.

pos = zeros(length(t),3);
for( ii=1:3 )
  pos(:,ii) = interp1(raystruct.time,raystruct.pos(:,ii),t,'pchip');
end;
vprel = zeros(length(t),3);
for( ii=1:3 )
  vprel(:,ii) = interp1(raystruct.time,raystruct.vprel(:,ii),t,'pchip');
end;
vgrel = zeros(length(t),3);
for( ii=1:3 )
  vgrel(:,ii) = interp1(raystruct.time,raystruct.vgrel(:,ii),t,'pchip');
end;
n = zeros(length(t),3);
for( ii=1:3 )
  n(:,ii) = interp1(raystruct.time,raystruct.n(:,ii),t,'pchip');
end;
B0 = zeros(length(t),3);
for( ii=1:3 )
  B0(:,ii) = interp1(raystruct.time,raystruct.B0(:,ii),t,'pchip');
end;
qs = zeros(length(t),raystruct.Nspec);
for( ii=1:raystruct.Nspec )
  qs(:,ii) = interp1(raystruct.time,raystruct.qs(:,ii),t,'pchip');
end;
ms = zeros(length(t),raystruct.Nspec);
for( ii=1:raystruct.Nspec )
  ms(:,ii) = interp1(raystruct.time,raystruct.ms(:,ii),t,'pchip');
end;
Ns = zeros(length(t),raystruct.Nspec);
for( ii=1:raystruct.Nspec )
  Ns(:,ii) = interp1(raystruct.time,raystruct.Ns(:,ii),t,'pchip');
end;
nus = zeros(length(t),raystruct.Nspec);
for( ii=1:raystruct.Nspec )
  nus(:,ii) = interp1(raystruct.time,raystruct.nus(:,ii),t,'pchip');
end;
