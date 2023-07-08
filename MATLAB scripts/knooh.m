% calculate rate of OH + NO -> HONO

% Updated 7/18/06 AEP

% Based on JPL Data Evaluation #15

%  rate=knooh(T,M)

function j=knooh(T,M)

k2o=7e-31.*(T./300).^(-2.6).*M;

k2oo=3.6e-11.*(T./300).^(-0.1);

k2=0.6.^((1+(log10(k2o./k2oo)).^2).^(-1.0));

l=real(k2);

j=(k2o./(1+(k2o./k2oo))).*l;