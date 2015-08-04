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
%          'sse2'   use SSE2 instructions
%          'avx'    use AVX  instructions
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
%   See also TETRACC, TOSYMMAT, SUB2UTM.
%
%   File    : pcc.m
%   Author  : Kristian Loewe

if nargin == 4                          % --- get user-specified settings
  tile = varargin{2};                   % tiling
  nthd = varargin{3};                   % number of threads
  impl = varargin{1};                   % implementation

elseif nargin == 1                      % --- auto-determine settings
  tile = 1;                             % tiling -> use COBL
  try                                   % nthd   -> try to determine #cores
    nthd = feature('numcores');
  catch ME
    nthd = round(nproc/2);
  end
  if     cpuinfo('avx')                 % impl.  -> depends on the cpu
    impl = 'avx';
  elseif cpuinfo('sse2')
    impl = 'sse2';
  else
    impl = 'naive';
  end

else
  error('Unexpected number of input arguments.');
end

variant = uint32(0);                    % --- check & parse settings
assert(isnumeric(tile) ...              % tiling
  && isscalar(tile) && tile >= 0);
tile = uint32(tile);
if     tile == 1                        % use COBL
  variant = bitor(variant, uint32(32));
elseif tile > 1                         % use standard tiling
  variant = bitor(variant, uint32(16));
end
assert(isnumeric(nthd) ...              % number of threads
  && isscalar(nthd) && nthd >= 0 && nthd <= 48);
nthd = uint32(nthd);
if nthd > 0
  variant = bitor(variant, uint32(64)); % enable multi-threading
end
assert(ischar(impl) ...                 % implementation
  && ismember(impl, {'avx','sse2','naive'}));
switch variant
  case 'avx'
    variant = bitor(variant, uint32(3));
  case 'sse2'
    variant = bitor(variant, uint32(2));
  case 'naive'
    variant = bitor(variant, uint32(1));
end

assert(ndims(data) == 2, ...            % --- check input data
  'Data array has unexpected number of dimensions.');
assert(isnumeric(data), ...
  'Data array is not numeric.');
assert(isreal(data), ...
  'Data array is not real.');

if strcmp(class(data),'single')         % --- call C/MEX program
  if tile == 1
    if nthd > 0                         % COBL tiling     / multiple threads
      r = mxPccXxFlt(data,variant,nthd);
    else                                % COBL tiling     / single thread
      r = mxPccXxFlt(data,variant);
    end
  elseif tile > 1
    if nthd > 0                         % standard tiling / multiple threads
      r = mxPccXxFlt(data,variant,tile,nthd);
    else                                % standard tiling / single thread
      r = mxPccXxFlt(data,variant,tile);
    end
  else
    if nthd > 0                         % no tiling       / multiple threads
      r = mxPccXxFlt(data,variant,nthd);
    else
      r = mxPccXxFlt(data,variant);
    end
  end
elseif strcmp(class(data),'double')
  if tile == 1
    if nthd > 0                         % COBL tiling     / multiple threads
      r = mxPccXxDbl(data,variant,nthd);
    else                                % COBL tiling     / single thread
      r = mxPccXxDbl(data,variant);
    end
  elseif tile > 1
    if nthd > 0                         % standard tiling / multiple threads
      r = mxPccXxDbl(data,variant,tile,nthd);
    else                                % standard tiling / single thread
      r = mxPccXxDbl(data,variant,tile);
    end
  else
    if nthd > 0                         % no tiling       / multiple threads
      r = mxPccXxDbl(data,variant,nthd);
    else
      r = mxPccXxDbl(data,variant);
    end
  end
else
  error('pcc:check_args','Data is of unexpected type.');
end

end
