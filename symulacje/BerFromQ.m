function [BER] = BerFromQ(Q)
%BERFROMQ Calculates BER from Quality factor
%
BER = 0.5*(1-erf(Q/sqrt(2)));

end

