function S = flip_phase(x)
% Samuel Wagner, UC Davis ECE MML, 2021
% flip_phase(x) returns a complex x with its
% angle negated, useful for hermetian construction

m = abs(x);                     % magnitude
a = -angle(x);                  % angle
S = m.*cos(a) + 1j.*m.*sin(a);  % flipped