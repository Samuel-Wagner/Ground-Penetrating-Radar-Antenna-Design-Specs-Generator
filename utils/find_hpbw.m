function hpbw = find_hpbw(patterndB,theta)
% Samuel Wagner, UC Davis ECE MML, 2021
% 
% function to calculate the half-power beamwidth of an antenna pattern
% pattern must be symmetric
% it will find the maximum index and traverse to the right
% after finding the half-power index, it calculates the angle and returns
% twice that value.

% inputs - patterndB - antenna pattern in dB
% theta  - angles at which patterndB is specified

[pattern_max,ind_max] = max(patterndB);			% maximum value & index
patterndB 	= patterndB-pattern_max;			% normalize the pattern
thresh		=	0.2;							% how close should it be to 3 dB?

% first index AFTER max at which the pattern is at least thresh near -3 dB
ind_right 	= find(abs(patterndB(ind_max:end) + 3) < thresh, 1,'first');	
ind_right 	= ind_max + ind_right - 1;	% re-shift to real value

% if find() returns an empty array, return at least a real number
if(isempty(ind_right))
    ind_right = numel(theta);
end

% assume pattern is symmetric, for simplicity.
hpbw = 2*theta(ind_right);