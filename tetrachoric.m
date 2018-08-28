function rt = tetrachoric(varargin)
%TETRACHORIC Estimate the tetrachoric correlation coefficient.
%
%   rt = tetrachoric(p11,p12,p21,p22) estimates the tetrachoric correlation
%   coefficient rt based on the specified cell probabilities p11,p12,p21,p22.
%   In the case of equal marginals, the formula rt = -cos(2*pi*p11) is used,
%   otherwise the 'bonett' algorithm (Bonett, 2006) is used, unless a
%   different method is specified (see below).
%
%   Using rt = tetrachoric(p) the cell probabilities are be specified by p,
%   where p is the four-element vector [p11,p12,p21,p22].
%
%   Using rt = tetrachoric(...,method) the user can specify which method
%   should be used to estimate rt. Available methods are 'cospi', 'edwards'
%   (Edwards and Edwards, 1984), and 'bonett' (Bonett, 2006). The default is
%   'bonett'.
%
%   Using rt = tetrachoric(...,'bonett',acc) the user can specify the
%   degree of accuracy to be used in Bonnet's algorithm (default: 0.0001).
%
%   Examples
%   --------
%   tetrachoric(1/3,1/6,1/6,1/3)               reports 0.5
%   tetrachoric(0.45,0.05,0.05,0.45)           reports 0.9511
%   tetrachoric(0.08,0.02,0.42,0.48)           reports 0.4395
%   tetrachoric(0.0967,0.0033,0.4033,0.4967)   reports 0.7492
%   tetrachoric(0.002,0.018,0.498,0.482)       reports -0.4734
%
%   p = [0.08,0.02,0.42,0.48];
%   tetrachoric(p, 'cospi')                    reports 0.5394
%   tetrachoric(p, 'edwards')                  reports 0.5348
%   tetrachoric(p, 'bonett')                   reports 0.4395
%
%   p = [0.3,0.1,0.2,0.4];
%   tetrachoric(p, 'cospi')                    reports 0.6132
%   tetrachoric(p, 'edwards')                  reports 0.6067
%   tetrachoric(p, 'bonett')                   reports 0.6071
%   format long g
%   tetrachoric(p, 'bonett', 0.1)              reports 0.59375
%   tetrachoric(p, 'bonett', 0.01)             reports 0.60546875
%   tetrachoric(p, 'bonett', 0.001)            reports 0.60693359375
%   tetrachoric(p, 'bonett', 0.0001)           reports 0.607086181640625
%   tetrachoric(p, 'bonett', 0.00001)          reports 0.607074737548828
%
%   References
%   ----------
%   Edwards, J.H. and Edwards, A.W.F.
%   Approximating the Tetrachoric Correlation Coefficient
%   Biometrics, Vol. 40, No. 2 (Jun., 1984), p. 563
%   https://www.jstor.org/stable/2531412
%
%   Bonett, D.G.
%   A new algorithm for the tetrachoric correlation
%   In I. J. Good (Section Editor), Comments, Conjectures and Conclusions
%   Journal of Statistical Computation and Simulation (2006) 76:8, 737-739.
%   https://dx.doi.org/10.1080/10629360500108186
%
%   Author: Kristian Loewe

method = 'bonett';
acc = 0.0001;

%%
assert(nargin >= 1 && nargin <= 6, 'Unexpected number of input arguments.');

if ~isscalar(varargin{1})   % rt = tetrachoric(p,...)
  assert(numel(varargin{1}) == 4);
  p11 = varargin{1}(1);
  p12 = varargin{1}(2);
  p21 = varargin{1}(3);
  p22 = varargin{1}(4);
else                        % rt = tetrachoric(p11,p12,p21,p22,...)
  p11 = varargin{1};
  p12 = varargin{2};
  p21 = varargin{3};
  p22 = varargin{4};
end

if nargin == 2 || nargin == 3
  assert(ischar(varargin{2}));
  method = varargin{2};
  if nargin == 3
    assert(isscalar(varargin{3}));
    acc = varargin{3};
  end

elseif nargin == 5 || nargin == 6
  assert(ischar(varargin{5}));
  method = varargin{5};
  if nargin == 6
    assert(isscalar(varargin{6}));
    acc = varargin{6};
  end
end

if (p11 + p12 + p21 + p22 > 1)
  n = p11 + p12 + p21 + p22;
  p11 = p11/n;
  p12 = p12/n;
  p21 = p21/n;
  p22 = p22/n;
end

assert(abs(p11 + p12 + p21 + p22 - 1) < 10*eps);

%%
if p11 == p22 && p12 == p21
  rt = -cos(2*pi*p11);

else
  switch method
    case 'cospi'
      rt = cos(pi/(1 + sqrt((p11*p22)/(p12*p21))));

    case 'edwards'
      a = ((p11*p22)/(p12*p21))^(pi/4);
      rt = (a - 1)/(a + 1);

    case 'bonett'
      pr = p11 + p12;
      pc = p11 + p21;
      h = norminv(pr);
      k = norminv(pc);
      L = -1;
      U = 1;
      er = 1;
      while (er > acc)
        rt = (L + U)/2;
        % p0 = mvncdf([h k], [0 0], [1 rt; rt 1]);
        % p0 = mvncdf([h k], [], [1 rt; rt 1]);
        p0 = internal.stats.bvncdf([h k], rt, 1e-8);
        p0(p0 < 0) = 0;
        p0(p0 > 1) = 1;
        er = U - L;
        if p0 > p11
          U = rt;
        else
          L = rt;
        end
      end
  end
end

end
