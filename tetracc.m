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
%      Implementation
%      v  'lut16'
%         'sse2'
%         'ssse3'
%         'pop32'
%         'pop64'
%         'm128i'
%
%      Tiling
%      t  0        none
%         2-n      standard tiling using t's next power of 2 as tile size
%
%      Multi-threading
%      p  0        auto  multi-threaded version (auto)
%                        the number of threads is automatically determined
%         1        none  single-threaded version
%         2-n      user  multi-threaded version (user)
%                        the number of threads is user-specified (through p)
%
%   Each variable in X is binarized with respect to its sample median; 
%   then for each pair of variables the tetrachoric correlation rt is 
%   found using rt = -cos(2*pi*n11/M), where n11 is the number of 
%   observations where both binarized variables are 1.
%
%   If you use this program, please cite:
%   Fast construction of voxel-level functional connectivity graphs
%   K. Loewe, M. Grueschow, C. Stoppel, R. Kruse, and C. Borgelt
%   BMC Neuroscience 15:78 (2014)
%
%   See also PCC, TOSYMMAT, SUB2UTM.
%
%   Filename : tetracc.m
%   Author   : Kristian Loewe

if nargin == 4                       % --- get user-specified settings
  impl = varargin{1};                % implementation (tetracc)
  tile = varargin{2};                % tiling
  nthd = varargin{3};                % multi-threading (number of threads)

elseif nargin == 1                   % --- auto-determine settings
  if hasIsaExtension('popcnt')       % implementation -> depends on the cpu
    impl = 'm128i';
  else
    impl = 'lut16';
  end
  tile = 0;                          % tiling -> none
  nthd = 0;                          % nthd   -> will be auto-determined

else
  error('Unexpected number of input arguments.');
end

variant = uint32(0);                 % --- check & parse settings
assert(ischar(impl) ...              % variant
  && ismember(impl, {'lut16','sse2','ssse3','pop32','pop64','m128i'}));
switch impl                          %   adjust variant code
  case 'lut16'
    variant = bitor(variant, uint32(1));
  case 'sse2'
    variant = bitor(variant, uint32(2));
  case 'ssse3'
    variant = bitor(variant, uint32(3));
  case 'pop32'
    variant = bitor(variant, uint32(4));
  case 'pop64'
    variant = bitor(variant, uint32(5));
  case 'm128i'
    variant = bitor(variant, uint32(6));
end
assert(isnumeric(tile) ...           % tiling
  && isscalar(tile) && tile >= 0 && tile ~= 1);
tile = uint32(tile);
if tile > 1                          %   enable standard tiling (variant code)
  variant = bitor(variant, uint32(16));
end
assert(isnumeric(nthd) ...           % multi-threading
  && isscalar(nthd));
if nthd == 0                         %   auto-determine number of threads
  try
    nthd = feature('numcores');
  catch ME %#ok<NASGU>
    nthd = round(nproc()/2);
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
  if tile > 1
    if nthd > 0                      % standard tiling / multiple threads
      r = mxTetraccXxFlt(data,variant,tile,nthd);
    else                             % standard tiling / single thread
      r = mxTetraccXxFlt(data,variant,tile);
    end
  else
    if nthd > 0                      % no tiling       / multiple threads
      r = mxTetraccXxFlt(data,variant,nthd);
    else                             % no tiling       / single thread
      r = mxTetraccXxFlt(data,variant);
    end
  end
elseif strcmp(class(data),'double')
  if tile
    if nthd > 0                      % standard tiling / multiple threads
      r = mxTetraccXxDbl(data,variant,tile,nthd);
    else                             % standard tiling / single thread
      r = mxTetraccXxDbl(data,variant,tile);
    end
  else
    if nthd > 0                      % no tiling       / multiple threads
      r = mxTetraccXxDbl(data,variant,nthd);
    else                             % no tiling       / single thread
      r = mxTetraccXxDbl(data,variant);
    end
  end
else
  error('Data array is of unexpected type.');
end

end

