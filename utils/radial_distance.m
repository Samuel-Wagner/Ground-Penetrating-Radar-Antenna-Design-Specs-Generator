function r = radial_distance(r1,r2,deassert)
    if(~exist('deassert','var'))
        deassert=false;
    end
    if(~deassert)
        %calculates radial distance between two points r1 and r2
        assert(numel(r1)==3,"r1 must have three elements");
        assert(numel(r2)==3,"r2 must have three elements");
        assert(ismember(string(class(r1)), ["double","single"]),"r1 must be of type double or single");
        assert(ismember(string(class(r2)), ["double","single"]),"r2 must be of type double or single");
    end
    r = sqrt((r1(1)-r2(1))^2+(r1(2)-r2(2))^2+(r1(3)-r2(3))^2);
end

