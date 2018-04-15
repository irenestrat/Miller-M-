function [h] = rayleigh_fading(p)

a_1 = abs(sqrt(0.5) * (randn + 1i*randn));
a_2 = abs(sqrt(0.5) * (randn + 1i*randn));

phi_1 = unifrnd(0, 2*pi);
phi_2 = unifrnd(0, 2*pi);

h_1 = a_1*exp(-1i*phi_1);
h_2 = a_2*exp(-1i*phi_2);

if (p == true)
    h = h_1*h_2;
else
    h = h_1;
end

end