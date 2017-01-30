function f = logsumexp(x)
% performs logsumexp calculation
%
% INPUT:
% x - vector of values
%
% OUTPUT:
% f - logsumexp calculation

the_max  = max(x);
x = x - the_max;
f = the_max + log(sum(exp(x)));

end