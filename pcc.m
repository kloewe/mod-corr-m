function r = pcc(data,varargin)
%PCC Compute pairwise Pearson correlation coefficients.
%   R = PCC(X) returns a vector of pairwise correlation coefficients
%   for an M-by-N data array X. Rows and columns of X correspond to
%   observations and variables, respectively. In other words, there
%   are M observations for each of the N variables in X.
%   The N*(N-1)/2 values in R correspond to the upper triangular
%   part of the symmetric N-by-N matrix Z of correlation coefficients.
%   Z can be obtained using Z = toSymMat(R). The correlation between
%   the Ith and Jth variable in X can be acessed using
%   R((I-1)*(N-I/2)+J-I) or Z(I,J).
%
%   Using R = PCC(X) the implementation variant to be used is
%   automatically chosen, depending on the features of your system.
%
%   Alternatively, the implementation variant to be used can be
%   specified by the user using R = PCC(X,v,t,p)
%
%      v   'naive'
%          'sse2'
%          'avx'
%
%      t   0        no tiling
%          1        cache-oblivious tiling
%          2-n      standard tiling using t's next power of 2 as tile size
%
%      p   0        no threads
%          2-n      parallelization using p threads
%
%   If you use this program, please cite:
%   Fast construction of voxel-level functional connectivity graphs
%   K. Loewe, M. Grueschow, C. Stoppel, R. Kruse, and C. Borgelt
%   BMC Neuroscience 15:78 (2014)
%
%   See also TETRACC, TOSYMMAT.
%
%   File    : pcc.m
%   Author  : Kristian Loewe

if nargin == 4                         % if optional arguments were passed

  variant = varargin{1};               % variant: get passed value
  assert(ischar(variant));             % make sure passed value is char
  vals = {'naive','sse2','avx'};       % define valid values
  assert(any(strcmp(variant,vals)));   % make sure passed value is valid

  tile = uint32(varargin{2});          % tiling:  get passed value
  assert(tile >= 0 ...                 % make sure passed value is 
    && (round(tile) == tile));         % a natural number >= 0

  nthreads = uint32(varargin{3});      % threads: get passed value
  assert(nthreads >= 0 ...             % make sure passed value is 
    && (round(nthreads) == nthreads)); % a natural number >= 0

  switch variant                       % select variant
    case 'avx'
      variant = uint32(3);             % avx
    case 'sse2'            
      variant = uint32(2);             % sse2
    case 'naive'
      variant = uint32(1);             % naive
  end

  if tile == 1                         % use cache-oblivious tiling
    variant = bitor(variant,uint32(32));
  elseif tile > 1                      % use standard tiling
    variant = bitor(variant,uint32(16));
  end

  if nthreads > 0                      % use threads
    variant = bitor(variant,uint32(64));
  end

elseif nargin == 1                     % determine implementation variant
  variant  = uint32(0);                % automatically (on the C side)
  tile     = uint32(0);
  nthreads = uint32(0);

else
  error('pcc:check_args', 'Unexpected number of input arguments.');
end

assert(isnumeric(data), ...            % make sure data is numeric ...
  'pcc:check_args', 'Data is not numeric.');
assert(isreal(data), ...               % ... and not complex
  'pcc:check_args', 'Data is not real.');

if strcmp(class(data),'single')        % call the appropriate mex/c - program
  if tile == 1
    if nthreads > 0
      r = mxPccXxFlt(data,variant,nthreads);      % cobl. tiling / threads
    else
      r = mxPccXxFlt(data,variant);               % cobl. tiling / no thrds
    end
  elseif tile > 1
    if nthreads > 0
      r = mxPccXxFlt(data,variant,tile,nthreads); % std. tiling  / threads
    else
      r = mxPccXxFlt(data,variant,tile);          % std. tiling  / no thrds
    end
  else
    if nthreads > 0
      r = mxPccXxFlt(data,variant,nthreads);      % no tiling    / threads
    else
      r = mxPccXxFlt(data,variant);
    end
  end
elseif strcmp(class(data),'double')
  if tile == 1
    if nthreads > 0
      r = mxPccXxDbl(data,variant,nthreads);      % cobl. tiling / threads
    else
      r = mxPccXxDbl(data,variant);               % cobl. tiling / no thrds
    end
  elseif tile > 1
    if nthreads > 0
      r = mxPccXxDbl(data,variant,tile,nthreads); % std. tiling  / threads
    else
      r = mxPccXxDbl(data,variant,tile);          % std. tiling  / no thrds
    end
  else
    if nthreads > 0
      r = mxPccXxDbl(data,variant,nthreads);      % no tiling    / threads
    else
      r = mxPccXxDbl(data,variant);
    end
  end
else
  error('pcc:check_args','Data is of unexpected type.');
end

end
