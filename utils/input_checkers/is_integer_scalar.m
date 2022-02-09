function tf = is_integer_scalar(x)
    tf =  numel(x)==1 && isfinite(x) && x==floor(x);    
end

