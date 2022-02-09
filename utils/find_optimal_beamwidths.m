function [hpbw_min, hpbw_max] = find_optimal_beamwidths(HPBWs,HRs)
    
    % sort to be ascending beamwidths, just in case.
    [HPBWs, sort_inds] = sort(HPBWs,'ascend');
    HRs = HRs(sort_inds);
    
    % find the minimum index and value
    [HR_at_min, min_ind] = min(HRs);
    
    % take an additional 10% for security
    HR_at_edges = 1.1 * HR_at_min;
        
    % split the HRs into left & right
    left  = HRs(1:min_ind);
    right = flip(HRs(min_ind:end));
    
    % tolerance for errors
    tolerance = HR_at_edges/100;
    
    % find the first time that HRs is within tolerance of hpbw_at_edges
    left_ind  = find(abs(left-HR_at_edges)<tolerance,1,'first');
    right_ind = find(abs(right-HR_at_edges)<tolerance,1,'first');

    if(isempty(right_ind))
       right_ind = 1; 
    end
    if(isempty(left_ind))
       left_ind = 1; 
    end
    
    % 
    right_ind = (numel(right)-right_ind) + min_ind;

    hpbw_min = HPBWs(left_ind);
    hpbw_max = HPBWs(right_ind);
    
end