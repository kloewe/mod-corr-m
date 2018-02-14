%TEST_PCC
%
%   Requires the following modules:
%
%     cpuinfo-m
%     corr-m
%
%   Author: Kristian Loewe

if ~hasIsaExtension('avx')
  warning('testPcc:NoAvxSupport', ...
    'AVX-based variants cannot be tested on this computer.\n');
end

verbose = 0;
quick = 0;


%% test 1
fprintf('\nTest 1 ...\n');

% single
a = [0.4430098 0.5680599 0.3433261 0.9175944 0.3769222 0.9770406];
b = [0.1130724 0.1368889 0.5478951 0.0083151 0.9292043 0.1646702];
c1 = pcc(single([a(:), b(:)]));
c2 = corrcoef(single([a(:), b(:)]));
c2 = c2(2);
fprintf('single ->  pcc: %f  corrcoef: %f  diff: %e\n', c1, c2, c2-c1);
assert( abs(c2-c1) < 10*eps('single') );

% double
a = [0.325733108917976 0.060120769323398 0.306466746642906 0.568924908300412];
b = [0.606927621835917 0.747400804922340 0.171137060673563 0.490894896500361];
c1 = pcc(double([a(:), b(:)]));
c2 = corrcoef([a(:), b(:)]);
c2 = c2(2);
fprintf('double ->  pcc: %f  corrcoef: %f  diff: %e\n', c1, c2, c2-c1);
assert( abs(c2-c1) < 10*eps('double') );

fprintf('[PASSED]\n\n');


%% test 2 (actually a series of tests)
fprintf('Test 2 ...\n');
if quick
  nV = 5000;
  nT = 200;
else
  nV = [2,3,4,7,15,17,501,1000,5000];
  nT = [4,6,8,9,15,16,17,32,50,64,100,128,175,200,256,300,400,500,512]; %,2,3
  % TODO: check why AVX fails for T < 4
end

nP = 0:proccnt();

for dtype = {'single','double'}
  for iV = 1:numel(nV)
    for iT = 1:numel(nT)
      if verbose; fprintf('\nV: %u  T: %u\n', nV(iV), nT(iT)); end

      data = rand(nT(iT), nV(iV), dtype{1});      % generate random data
      limmad = 5*eps(dtype{1});
      limsad = limmad*(nV(iV)*(nV(iV)-1)/2);

      res1 = corrcoef(data);                      % corrcoef (!REFERENCE!)
      res1 = res1(logical(tril(ones(nV(iV)),-1)));
      res1 = res1(:);

      % --- single-thread variants ---
      if verbose; fprintf('naive\n'); end
      res2 = pcc(data, 'naive', 0, 0);            % naive
      mad = max(abs(res1-res2));                  % max abs diffs
      sad = sum(abs(res1-res2));                  % sum of abs diffs
      if verbose; fprintf('  mad:  %e\n', mad); end
      if verbose; fprintf('  sad:  %e\n', sad); end
      assert( mad < limmad && sad < limsad );
      res2b = pcc(data, 'naive', 2, 0);           % naive + tiling
      assert(isequal(res2,res2b));
      res2c = pcc(data, 'naive', 1, 0);           % naive + cobl
      assert(isequal(res2,res2c));

      if verbose; fprintf('sse2\n'); end
      res3 = pcc(data, 'sse2', 0, 0);             % sse2
      mad = max(abs(res2-res3));                  % max abs diffs
      sad = sum(abs(res2-res3));                  % sum of abs diffs
      if verbose; fprintf('  mad:  %e\n', mad); end
      if verbose; fprintf('  sad:  %e\n', sad); end
      assert( mad < limmad && sad < limsad );
      res3b = pcc(data, 'sse2', 2, 0);            % sse2  + tiling
      assert(isequal(res3,res3b));
      res3c = pcc(data, 'sse2', 1, 0);            % sse2  + cobl
      assert(isequal(res3,res3c));

      if hasIsaExtension('avx')
        if verbose; fprintf('avx\n'); end
        res4 = pcc(data, 'avx', 0, 0);            % avx
        mad = max(abs(res2-res4));                % max abs diffs
        sad = sum(abs(res2-res4));                % sum of abs diffs
        if verbose; fprintf('  mad:  %e\n', mad); end
        if verbose; fprintf('  sad:  %e\n', sad); end
        assert( mad < limmad && sad < limsad );
        res4b = pcc(data, 'avx', 2, 0);           % avx  + tiling
        assert(isequal(res4,res4b));
        res4c = pcc(data, 'avx', 1, 0);           % avx  + cobl
        assert(isequal(res4,res4c));
      end

      % --- multi-thread variants ---
      if ~quick
        for iP = 1:numel(nP)
          res2d = pcc(data, 'naive', 0, nP(iP));  % naive + threads
          assert(isequal(res2,res2d));
          res2e = pcc(data, 'naive', 2, nP(iP));  % naive + tiling + threads
          assert(isequal(res2,res2e));
          res2f = pcc(data, 'naive', 1, nP(iP));  % naive + cobl   + threads
          assert(isequal(res2,res2f));

          res3d = pcc(data, 'sse2',  0, nP(iP));  % sse2 + threads
          assert(isequal(res3,res3d));
          res3e = pcc(data, 'sse2',  2, nP(iP));  % sse2 + tiling  + threads
          assert(isequal(res3,res3e));
          res3f = pcc(data, 'sse2',  1, nP(iP));  % sse2 + cobl    + threads
          assert(isequal(res3,res3f));

          if hasIsaExtension('avx')
            res4d = pcc(data, 'avx',  0, nP(iP)); % avx + threads
            assert(isequal(res4,res4d));
            res4e = pcc(data, 'avx',  2, nP(iP)); % avx + tiling   + threads
            assert(isequal(res4,res4e));
            res4f = pcc(data, 'avx',  1, nP(iP)); % avx + cobl     + threads
            assert(isequal(res4,res4f));
          end
        end
      end

    end
  end
end
fprintf('[PASSED]\n\n');

%% test 3
fprintf('Test 3 ...\n');

for dtype = {'single','double'}
  for i = 1:10000;
    d = rand(200, 1, dtype{1});
    assert(pcc([d -d]) >= -1.0 & pcc([d d]) <= 1.0);
  end
end

fprintf('[PASSED]\n\n');
