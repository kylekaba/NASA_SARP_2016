% calculate rate of HO2 + HO2 -> H2O2 + O2
% rate=kho2ho2(T,M)
% uncertainties   A  1.3  E/R  +-200
%  p dependent    A  1.3   E/R  +-400
%IMP 8/4/05 update JPL 2003
function j=kho2ho2(T,M)
% j=1.7e-33.*M.*exp(1000./T);

R = 8.3144621;
j = 2.2*10^-13.*exp(4989./(R.*T));

% m = 1 + (1.40D-21*EXP(2200/TEMP)*H2O)
% 2.20D-13*m*EXP(600/TEMP)