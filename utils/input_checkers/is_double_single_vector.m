function tf = is_double_single_vector(x)
    tf =  (...
            strcmp(string(class(x)),"double") ...
            || strcmp(string(class(x)),"single")...
           )...
          && (numel(x)>1) && numel(size(x))==2 && abs(diff(size(x)))>0;
end