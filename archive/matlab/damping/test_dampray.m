const;
addpath('../../bin');
addpath('../');

TOL = 1e-3;
INTEGMETHOD = 'accurate';

f = 4000;
w = 2*pi*f;

% ray index
ri = 1;

out=readrayoutput('test_4000Hz.txt');
t = linspace(out{ri}.time(1),out{ri}.time(end),300).';
t = out{ri}.time(1:10:end);
%t = out{ri}.time;

[pos,vprel,vgrel,n,B0,qs,ms,Ns,nus] = interpolaterayoutput(out{ri},t);

% Ignore collisions for this code
nus = 0*nus;

% SET UP THE HOT ELECTRON DISTRIBUTION
% electron mass
me = ME;
% electron charge
qe = -Q;
% maxwellian distribution
% Temperature in eV 
%kT = 1e3*Q; % 1 keV
%kT = 2e3*Q; % 2 keV
%Ne_h = .1e6;
%fe = @(vperp,vpar) Ne_h*maxwellboltzmann(vperp,vpar,me,kT);
% Suprathermal distribution from Bell 2002
fe = @(vperp,vpar) suprathermal(vperp,vpar);
% Resonances (0=landau, +-1 = cyclotron)
m = [-1 0 1];
% cell array of hot distribution functions
fs = {fe};
  
magnitude = zeros(size(t));
magnitude(1) = 1;
kis = zeros(size(t));
disp('Computing damping coefficient along ray path...');
for( ii=2:length(t) )
  fprintf('%d of %d\n', ii-1, length(t)-1)
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
    % Take the component of ki along vg (yes, this is right!).
    ki_along_vg = ki*(k*vgrel(ii,:)')/(norm(k)*norm(vgrel(ii,:)));

    dist = norm(pos(ii,:)-pos(ii-1,:));
    kis(ii) = ki_along_vg;
    magnitude(ii) = magnitude(ii-1)*exp(-dist*ki_along_vg);
  else
    disp('Re{n} = 0, not solving evanescent mode');
  end;
end;
% Power
disp('Done');

cmap = jet(64);
minmag = min(magnitude);
maxmag = max(magnitude);
figure(1)
for( ii=2:length(t) )
  ind = floor(63*(magnitude(ii)-minmag)/(maxmag-minmag))+1;
  plot3(pos((ii-1):ii,1),pos((ii-1):ii,2),pos((ii-1):ii,3),...
        'Color', cmap(ind,:),'LineWidth',2);
  hold on
end;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
ax=axis;
pad = .025;
axis([ax(1)-pad*(ax(2)-ax(1)),...
      ax(2)+pad*(ax(2)-ax(1)),...
      ax(3)-pad*(ax(4)-ax(3)),...
      ax(4)+pad*(ax(4)-ax(3)),...
      ax(5)-pad*(ax(6)-ax(5)),...
      ax(6)+pad*(ax(6)-ax(5))...
     ]);
view(3)
h=colorbar;
set(h,'YTick',linspace(1,65,6));
len = length(get(h,'YTick'));
numbersnew = linspace(minmag,maxmag,len);
numberslabels = {};
for( ii=1:len )
  numberslabels{ii} = sprintf('%2.2g', numbersnew(ii));
end;
set(h,'YTickLabel',numberslabels);
title('Ray amplitude');


minmag = min(kis);
maxmag = max(kis);
figure(2)
for( ii=2:length(t) )
  ind = floor(63*(kis(ii)-minmag)/(maxmag-minmag))+1;
  plot3(pos((ii-1):ii,1),pos((ii-1):ii,2),pos((ii-1):ii,3),...
        'Color', cmap(ind,:),'LineWidth',2);
  hold on
end;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
ax=axis;
pad = .025;
axis([ax(1)-pad*(ax(2)-ax(1)),...
      ax(2)+pad*(ax(2)-ax(1)),...
      ax(3)-pad*(ax(4)-ax(3)),...
      ax(4)+pad*(ax(4)-ax(3)),...
      ax(5)-pad*(ax(6)-ax(5)),...
      ax(6)+pad*(ax(6)-ax(5))...
     ]);
view(3)
h=colorbar;
set(h,'YTick',linspace(1,65,6));
len = length(get(h,'YTick'));
numbersnew = linspace(minmag,maxmag,len);
numberslabels = {};
for( ii=1:len )
  numberslabels{ii} = sprintf('%2.2g', numbersnew(ii));
end;
set(h,'YTickLabel',numberslabels);
title('Damping ratio along v_g (m^{-1})');

ti=linspace(t(1), t(end),100 );
magnitude_i = interp1(t,magnitude,ti,'pchip');
figure(3)
hold on
plot(ti,magnitude_i.^2,'r');
title('Power');
xlabel('time (s)');
