% calculate rate of OH + NO2 + M -> HNO3 + M
%  rate=kohno2(T,M)
% SP 12/10 update Mollner et al. 2010, Science

function j=kohno2new(T,M)
k9o=1.48e-30.*(T./300).^(-3).*M;
k9oo=2.58e-11.*(T./300).^(0);
k9=(k9o./(1+(k9o./k9oo))).*0.6.^((1+(log10(k9o./k9oo)).^2).^(-1.0));
j=k9;
