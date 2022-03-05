function [Q] = CalculateQ(v1a,v1b,v0a,v0b)
%%%v1a
%%%v1b
%%%
%%%v0a
%%%v0b

s1 = (v1a-v1b)/2;
v1 = (v1a+v1b)/2;

s0 = (v0a-v0b)/2;
v0 = (v0a+v0b)/2;

Q=(v1-v0)/(s1+s0);
end

