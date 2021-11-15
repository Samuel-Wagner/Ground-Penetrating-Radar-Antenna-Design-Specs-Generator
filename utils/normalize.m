function n = normalize(N)
	% normalize
	% function to normalize an input N-D matrix N by its
	% largest absolute value
	%
	% N - input (un-normalized) matrix
	% n - output (normalized) matrix
	n = N ./ max(abs(N(:)));
end