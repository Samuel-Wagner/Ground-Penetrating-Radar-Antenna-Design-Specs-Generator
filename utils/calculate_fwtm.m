function fwtm = calculate_fwtm(t,x)
% Samuel Wagner, UC Davis ECE MML, 2021

% function to calculate the full-width at tenth-max of a pulse
% for simplicity, this function is designed for symmetric pulses
% won't work with pulses with very non-smooth enevelopes

% inputs
% t - time vector at which x is sampled (s)
% x - pulse signal vector

fwtm_threshold      = 0.1;              % when should the fwtm be "called"?
env                 = abs(hilbert(x));  % calculate the envelope of x
[max_mag,max_ind]   = max(env);         % index, value of max of envelope

right               = env(max_ind:end);   % right side of signal
ind_lt_tenth        = find(right/max_mag < fwtm_threshold, 1, 'first'); % when env falls below threshold

% return fwtm = 2 * distance from center to right edge (symmetric assumption)
fwtm = (t(ind_lt_tenth + max_ind - 1) - t(max_ind)) * 2; 