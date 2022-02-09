function fwtm = calculate_fwtm(t,x)
% function to calculate the full-width at
% tenth-max of a pulse
% this function is designed for symmetric pulses
% but might work with non-symmetric pulses if they fall off
% to the right side. won't work with pulses with very non-smooth enevelopes

% inputs
% t - time vector at which x is sampled
% x - signal, real

env                 = abs(hilbert(x));  % find the envelope
[max_mag,max_ind]   = max(env);         % index, value of max

right               = env(max_ind:end);   % right side of signal
ind_lt_tenth        = find(right/max_mag < 0.1, 1, 'first');

fwtm = (t(ind_lt_tenth + max_ind - 1) - t(max_ind)) * 2; 


