function k = kho2ro2(T,M)

%k=3.5e-13*exp(1000./T); % Taken from IUPAC/JPC Reference data
								% generic rate constant for ro2 + ho2 --> rooh + o2
                                
k = 2.9e-13.*exp(1300./T); % From MCM v3.1, 2003. -BWL
% R = 8.3144621;
% k = 3.8*10^-13*exp(6485./(R.*T));