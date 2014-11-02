function r = tetracc(data,varargin)
%TETRACC Compute pairwise tetrachoric correlation coefficients.
%   R = TETRACC(X) returns a vector of pairwise correlation coefficients
%   for an M-by-N data array X. Rows and columns of X correspond to
%   observations and variables, respectively. In other words, there
%   are M observations for each of the N variables in X.
%   The N*(N-1)/2 values in R correspond to the upper triangular
%   part of the symmetric N-by-N matrix Z of correlation coefficients.
%   Z can be obtained using using Z = toSymMat(R). The correlation between
%   the Ith and Jth variable in X can be acessed using
%   R((I-1)*(N-I/2)+J-I) or Z(I,J).
%
%   Using R = TETRACC(X) the implementation variant to be used is
%   automatically chosen, depending on the features of your system.
%
%   Alternatively, the implementation variant to be used can be
%   specified by the user using R = TETRACC(X,v,t,p)
%
%      v   'lut16'
%          'sse2'
%          'ssse3'
%          'pop32'
%          'pop64'
%          'm128i'
%
%      t   0        no tiling
%          2-n      standard tiling using t's next power of 2 as tile size
%
%      p   0        no threads
%          2-n      parallelization using p threads
%
%   Each variable in X is binarized with respect to its sample median; 
%   then for each pair of variables the tetrachoric correlation rt is 
%   found using rt = -cos(2*pi*n11/M), where n11 is the number of 
%   observations where both binarized variables are 1.
%
%   See also PCC, TOSYMMAT.
%
%   Filename : tetracc.m
%   Author   : Kristian Loewe

if nargin == 4                         % if optional arguments were passed

  variant = varargin{1};               % variant: get passed value
  assert(ischar(variant));             % make sure passed value is char
  vals = {'lut16','sse2', ...          % define valid values
    'ssse3','pop32', ...
    'pop64','m128i'};
  assert(any(strcmp(variant,vals)));   % make sure passed value is valid

  tile = uint32(varargin{2});          % tiling:  get passed value
  assert(tile >= 0 ...                 % make sure passed value is 
    && (round(tile) == tile));         % a natural number >= 0

  nthreads = uint32(varargin{3});      % threads: get passed value
  assert(nthreads >= 0 ...             % make sure passed value is 
    && (round(nthreads) == nthreads)); % a natural number >= 0

  switch variant                       % select variant
    case 'lut16'
      variant = uint32(1);             % lut16
    case 'sse2'
      variant = uint32(2);             % sse2
    case 'ssse3'
      variant = uint32(3);             % ssse3
    case 'pop32'
      variant = uint32(4);             % pop32
    case 'pop64'
      variant = uint32(5);             % pop64
    case 'm128i'
      variant = uint32(6);             % m128i
  end

  if tile > 0                          % use standard tiling
    variant = bitor(variant,uint32(16));
  end

  if nthreads > 0                      % use threads
    variant = bitor(variant,uint32(32));
  end

elseif nargin == 1                     % determine implementation variant
  variant  = uint32(0);                % automatically (on the C side)
  tile     = uint32(0);
  nthreads = uint32(0);

else
  error('tetracc:check_args', 'Unexpected number of input arguments.');
end

assert(isnumeric(data), ...            % make sure data is numeric ...
  'tetracc:check_args', 'Data is not numeric.');
assert(isreal(data), ...               % ... and not complex
  'tetracc:check_args', 'Data is not real.');

if strcmp(class(data),'single')        % call the appropriate mex/c - program
  if tile                                              % tiling threads
    if nthreads > 0
      r = mxTetraccXxFlt(data,variant,tile,nthreads);  %    1      1
    else
      r = mxTetraccXxFlt(data,variant,tile);           %    1      0
    end
  else
    if nthreads > 0
      r = mxTetraccXxFlt(data,variant,nthreads);       %    0      1
    else
      r = mxTetraccXxFlt(data,variant);
    end
  end
elseif strcmp(class(data),'double')
  if tile
    if nthreads > 0
      r = mxTetraccXxDbl(data,variant,tile,nthreads);  %    1      1
    else
      r = mxTetraccXxDbl(data,variant,tile);           %    1      0
    end
  else
    if nthreads > 0
      r = mxTetraccXxDbl(data,variant,nthreads);       %    0      1
    else
      r = mxTetraccXxDbl(data,variant);
    end
  end
else
  error('tetracc:check_args','Data is of unexpected type.');
end

end
