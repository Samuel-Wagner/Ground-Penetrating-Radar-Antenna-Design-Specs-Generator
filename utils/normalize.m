function n = normalize(N)
% Samuel Wagner, UC Davis ECE MML, 2021
% divides a matrix by its maximum magnitude element
n = N ./ max(abs(N(:)));