function j=kro2no(T,M)

R = 8.3144621;

% j=2.8e-12.*exp(2370../(R.*T));
 j=2.8e-12.*exp(285./T);
% j=2.8e-12.*exp(300./T);
% 2.8x10-12 [cm3/molecule s] e2370 [±831 J/mole]/RT
