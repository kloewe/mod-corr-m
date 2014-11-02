function matrix = toSymMat(vect,varargin)
%TOSYMMAT Format correlation matrix.
%   S = TOSYMMAT(R) converts R into a square, symmetric matrix S,
%   where R is a vector of pairwise correlation coefficients as
%   created by the functions PCC or TETRACC.
%   By default, S has ones along its diagonal.
%   Using S = TOSYMMAT(R,d), the diagonal values are replaced by d.
%
%   See also PCC, TETRACC.
%
%   File    : toSymMat.m
%   Author  : Kristian Loewe

vd = 1;                              % default val for diagonal elements

if nargin == 2                       % user-specified value for
  vd = varargin{1};                  % diagonal elements
end

vect = vect(:);

funcX = @(n) n*(n-1)/2;              % the number of entries x in the
                                     % upper/lower triangular part of an
                                     % n-by-n matrix is given by x = n(n-1)/2

funcN = @(x) (sqrt(8*x+1)+1)/2;      % n can be found given x using
                                     % n = (sqrt(8*x+1)+1)/2

x = numel(vect);
n = funcN(x);

assert(funcX(n) == x, ...
  'toSymMat', ...
  'Inconsistent sizes of input vector and output matrix');
assert(ceil(n) == floor(n), ...      % is n integer?
  'toSymMat', ...
  'Unexpected result when calculating size of output matrix.');


if isnumeric(vect)
  matrix = zeros(n,n,class(vect));
  matrix(tril(true(n,n),-1)) = vect;
  matrixT = matrix';
  matrix = matrix + matrixT;
  matrix(1:(n+1):n*n) = vd;          % set main diagonal elements to vd
  
elseif islogical(vect)
  matrix = false(n,n);
  matrix(tril(true(n,n),-1)) = vect;
  matrixT = matrix';
  matrix = matrix | matrixT;
  matrix(1:(n+1):n*n) = logical(vd); % set main diagonal elements to vd
end

end
