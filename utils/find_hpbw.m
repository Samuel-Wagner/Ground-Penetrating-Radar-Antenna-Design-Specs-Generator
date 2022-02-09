
function hpbw = find_hpbw(patterndB,theta)
    [~,ind_max] = max(patterndB);
    patterndB = patterndB-max(patterndB);
    thresh=0.2;
    ind_right = find(abs(patterndB(ind_max:end) +3 ) < thresh, 1,'first');
    ind_right = ind_max + ind_right - 1;
    
    if(isempty(ind_right))
        ind_right = numel(theta);
    end
    
    hpbw = 2*theta(ind_right);
end