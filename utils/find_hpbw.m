function hpbw = find_hpbw(patterndB,theta)
	% find_hpbw
	% function to find the half-power beamwidth of a given pattern
	%
	% inputs
	% patterndB - antenna gain pattern in dB
	% theta     - theta values at which patterndB is measured
	%
	% outputs
	% hpbw      - half-power beamwidth of patterndB

	thresh 		= 0.2;	% threshold to calculate HPBW lims
    [~,ind_max] = max(patterndB);
    patterndB 	= patterndB - max(patterndB);
    ind_right 	= find(abs(patterndB(ind_max:end) +3 ) < thresh, 1,'first');
    ind_right 	= ind_max + ind_right - 1;
    
    % just in case there is no HPBW.
    if(isempty(ind_right))
        ind_right = numel(theta);
    end
    
    hpbw = 2*theta(ind_right);
end