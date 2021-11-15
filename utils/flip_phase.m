function S = flip_phase(x)
    % flip_phase(x) returns a complex x with its
    % angle negated, useful for hermetian construction
    m = abs(x);                     % magnitude
    a = -angle(x);                  % angle
    S = m.*cos(a) + 1j.*m.*sin(a);  % flipped
end