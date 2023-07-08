% calculate rate of HO2 +NO -> OH + NO2
%  rate=kho2no(T,M)
%   uncertainties   A  1.2  E/R  +-80
% IMP 8/4/05 update JPL 2003
function j=kho2no(T,M)

R = 8.3144621;

 j=3.6e-12.*exp(2245../(R.*T));
% j=3.5e-12.*exp(250../(T));
% 3.6x10-12 [cm3/molecule s] e2245 [±831 J/mole]/RT
