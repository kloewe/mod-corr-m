%TESTTETRACC
%
%   Requires the following modules:
%
%     cpuinfo-m
%     corr-m
%
%   Author: Kristian Loewe

if ~hasIsaExtension('popcnt')
    warning('testTetracc:NoPopcntSupport', ...
        'POPCNT-based variants cannot be tested on this computer.\n');
end

quick = 0;


%% test 1
fprintf('Test 1 ...\n');
a = logical([1 1 1 0 0 0 1 0]);
b = logical([0 1 1 0 0 0 1 1]);
n11 = sum(a & b);
rt = -cos(2*pi*n11/8);
rt = single(rt);

aflt = single([5 5 5 2 2 2 5 2]);
bflt = single([2 5 5 2 2 2 5 5]);
% dichotomize([5 5 5 2 2 2 5 2]) -> [1 1 1 0 0 0 1 0]
% dichotomize([2 5 5 2 2 2 5 5]) -> [0 1 1 0 0 0 1 1]

assert(isequal(rt, tetracc([aflt(:), bflt(:)])));

fprintf('[PASSED]\n');


%% test 2 (actually a series of tests)
fprintf('Test 2 ...\n');
if quick
  nV = 5000;
  nT = 200;
else
  nV = [2,3,4,7,15,17,501,1000,5000];
  nT = [2:50,64,100,128,175,200,256,300,400,500,512];
end

nP = 0:proccnt();

for iV = 1:numel(nV)
  for iT = 1:numel(nT)
    data = rand(nT(iT),nV(iV),'single');    % generate random data

    res1 = tetracc(data, 'lut16', 0, 0);    % lut16 (!REFERENCE!)

    for tile = [0, 2.^[1:14]]
      for iP = 1:numel(nP)
        res2 = tetracc(data, 'lut16', tile, nP(iP));     % lut16
        assert(isequal(res1,res2));
        res2 = tetracc(data, 'sse2',  tile, nP(iP));     % sse2
        assert(isequal(res1,res2));
        res2 = tetracc(data, 'ssse3', tile, nP(iP));     % ssse3
        assert(isequal(res1,res2));
        if hasIsaExtension('popcnt')
          res2 = tetracc(data, 'pop32', tile, nP(iP));   % pop32
          assert(isequal(res1,res2));
          res2 = tetracc(data, 'pop64', tile, nP(iP));   % pop64
          assert(isequal(res1,res2));
          res2 = tetracc(data, 'm128i', tile, nP(iP));   % m128i
          assert(isequal(res1,res2));
        end
      end
    end
  end
end
fprintf('[PASSED]\n');

%% test 3:  Check against old version 
% nV = [2,501,1000,5000];
% nT = [3,50,100,175,200,300,400,500,2,4,8,16,32,64,128,256,512];
% for i = 1:10000
% for iV = 1:numel(nV)
%   for iT = 1:numel(nT)
%     data = rand(nT(iT),nV(iV),'single');
%     
%     addpath /opt/klt/fc-tools/fc-tools-0.1-38/simimat/
%     addpath /opt/klt/fc-tools/fc-tools-0.1-38/util/
%     addpath /opt/klt/fc-tools/fc-tools-0.1-38/dicho/
% %     which tetracc
%     res1 = tetracc(data,'parallel',0);
%     
%     
%     addpath /home/kris/workspace/projects/klt/fc-tools/fc-tools/m/corr/tbx/
% %     which tetracc
%     res2 = tetracc(data);
%     
% %     assert(isequal(res1,res2));
%     mad = max(abs(res1-res2));              % max abs diff
%     sad = sum(abs(res1-res2));              % sum of abs diffs
% %     if verbose
%       fprintf('\n  max of abs differences:  %e\n', mad);
%       fprintf('  sum of abs differences:  %e\n', sad);
% %     end
%     assert( mad < 1e-06 && sad < 5e-07*(nV(iV)*(nV(iV)-1)/2) );
%   end
% end
% end
% fprintf('[PASSED]\n');
% 

