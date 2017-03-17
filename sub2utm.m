function idx = sub2utm(I,J,n,first)
%SUB2UTM Linear index (wrt upper triangle) from matrix subscripts.
%   SUB2UTM is used to determine the equivalent linear index (wrt a
%   vector of upper triangular matrix elements) corresponding to the
%   specified row and column subscripts (wrt a symmetric n-by-n matrix).
%
%   idx = sub2utm(I,J,n)
%   idx = sub2utm(I,J,n,first)
%
%   first   1 (default)  1-based indexing (of elements and rows/cols)
%           0            0-based indexing (of elements and rows/cols)
%
%   See also TOSYMMAT, PCC, TETRACC.
%
%   File    : sub2utm.m
%   Author  : Kristian Loewe

assert(isscalar(n));
assert(isvector(I) && isvector(J));
assert(isequal(size(I),size(J)));
assert(~any(isnan(I)) && ~any(isnan(J)));
assert(~any(I==J));

if ~exist('first','var')
  first = 1;
end
assert(ismember(first,[0 1]));

if first == 1
  I = I-1;
  J = J-1;
end
assert( min(I) >= 0 && max(I) <= n-1 && min(J) >= 0 && max(J) <= n-1 );

idx = zeros(size(I));
for k = 1:numel(I)
  if I(k) > J(k)
    tmp = I(k);
    I(k) = J(k);
    J(k) = tmp;
  end
  idx(k) = I(k).*(n + n - I(k) - 3) ./ 2 - 1 + J(k);
end

if first == 1
  idx = idx + 1;
end

end
