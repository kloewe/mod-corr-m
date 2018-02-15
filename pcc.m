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
%      Vectorization
%      v  'naive'  none
%         'sse2'   use SSE2 instructions
%         'avx'    use AVX  instructions
%
%      Tiling
%      t  0        none
%         1        cache-oblivious tiling
%         2-n      standard tiling using t's next power of 2 as tile size
%
%      Multi-threading
%      p  0        auto  multi-threaded version (auto)
%                        the number of threads is automatically determined
%         1        none  single-threaded version
%         2-n      user  multi-threaded version (user)
%                        the number of threads is user-specified (through p)
%
%   If you use this program, please cite:
%   Fast construction of voxel-level functional connectivity graphs
%   K. Loewe, M. Grueschow, C. Stoppel, R. Kruse, and C. Borgelt
%   BMC Neuroscience 15:78 (2014)
%
%   See also TETRACC, TOSYMMAT, SUB2UTM.
%
%   File    : pcc.m
%   Author  : Kristian Loewe

if nargin == 4                       % --- get user-specified settings
  vect = varargin{1};                % vectorization
  tile = varargin{2};                % tiling
  nthd = varargin{3};                % multi-threading (number of threads)

elseif nargin == 1                   % --- auto-determine settings
  if     hasIsaExtension('avx')      % vect.  -> depends on the cpu
    vect = 'avx';
  elseif hasIsaExtension('sse2')
    vect = 'sse2';
  else
    vect = 'naive';
  end
  tile = 1;                          % tiling -> use COBL
  nthd = 0;                          % nthd   -> will be auto-determined

else
  error('Unexpected number of input arguments.');
end

variant = uint32(0);                 % --- check & parse settings
assert(ischar(vect) ...              % vectorization
  && ismember(vect, {'avx','sse2','naive'}));
if size(data,1) < 8                  % if M < 8 fallback to SSE2
  vect = 'sse2';                     % (todo: fix problem with AVX-version)
end
switch vect                          % adjust variant code
  case 'naive'
    variant = bitor(variant, uint32(1));
  case 'sse2'
    variant = bitor(variant, uint32(2));
  case 'avx'
    variant = bitor(variant, uint32(3));
end
assert(isnumeric(tile) ...           % tiling
  && isscalar(tile) && tile >= 0);
tile = uint32(tile);
if     tile == 1                     %   enable COBL (variant code)
  variant = bitor(variant, uint32(32));
elseif tile > 1                      %   enable standard tiling (variant code)
  variant = bitor(variant, uint32(16));
end
assert(isnumeric(nthd) ...           % multi-threading
  && isscalar(nthd));
if nthd == 0                         %   auto-determine number of threads
  ncores = corecnt();
  nprocs = proccnt();
  if ncores == nprocs/2
    nthd = ncores + ncores/2;
  else
    nthd = ncores;
  end
  assert(isnumeric(nthd) && isscalar(nthd) && nthd >= 0 && nthd <= 1024);
end
nthd = uint32(nthd);
if nthd > 0                          %   enable multi-threading (variant code)
  variant = bitor(variant, uint32(64));
end

assert(ndims(data) == 2, ...         % --- check input data
  'Data array has unexpected number of dimensions.');
assert(isnumeric(data), ...
  'Data array is not numeric.');
assert(isreal(data), ...
  'Data array is not real.');
assert(~isempty(data), ...
  'Data array is empty.');

if strcmp(class(data),'single')      % --- call the appropriate C/MEX routine
  if tile == 1
    if nthd > 0                      % COBL tiling     / multiple threads
      r = mxPccXxFlt(data,variant,nthd);
    else                             % COBL tiling     / single thread
      r = mxPccXxFlt(data,variant);
    end
  elseif tile > 1
    if nthd > 0                      % standard tiling / multiple threads
      r = mxPccXxFlt(data,variant,tile,nthd);
    else                             % standard tiling / single thread
      r = mxPccXxFlt(data,variant,tile);
    end
  else
    if nthd > 0                      % no tiling       / multiple threads
      r = mxPccXxFlt(data,variant,nthd);
    else                             % no tiling       / single thread
      r = mxPccXxFlt(data,variant);
    end
  end
elseif strcmp(class(data),'double')
  if tile == 1
    if nthd > 0                      % COBL tiling     / multiple threads
      r = mxPccXxDbl(data,variant,nthd);
    else                             % COBL tiling     / single thread
      r = mxPccXxDbl(data,variant);
    end
  elseif tile > 1
    if nthd > 0                      % standard tiling / multiple threads
      r = mxPccXxDbl(data,variant,tile,nthd);
    else                             % standard tiling / single thread
      r = mxPccXxDbl(data,variant,tile);
    end
  else
    if nthd > 0                      % no tiling       / multiple threads
      r = mxPccXxDbl(data,variant,nthd);
    else                             % no tiling       / single thread
      r = mxPccXxDbl(data,variant);
    end
  end
else
  error('Data array is of unexpected type.');
end

end

