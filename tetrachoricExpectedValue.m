function Er = tetrachoricExpectedValue(p11,p12,p21,p22,N,useConstraints)
%TETRACHORICEXPECTEDVALUE Expected value of the tetrachoric correlation.
%
%   Given the cell probabilities p11, p12, p21, p22 and the total frequency N,
%   Er = tetrachoricExpectedValue(p11,p12,p21,p22,N) computes the expected
%   value of the tetrachoric correlation coefficient based on equations 3
%   and 4 in the paper by Brown and Benedetti (1977).
%
%   Limitations
%   -----------
%   N needs to be even.
%
%   Examples
%   ---------
%   tetrachoricExpectedValue(1/3,1/6,1/6,1/3,30)
%   tetrachoricExpectedValue(1/3,1/6,1/6,1/3,90)
%
%   tetrachoricExpectedValue(0.45,0.05,0.05,0.45,40)
%   tetrachoricExpectedValue(0.45,0.05,0.05,0.45,100)
%
%   tetrachoricExpectedValue(0.08,0.02,0.42,0.48,100)
%
%   tetrachoricExpectedValue(0.075,0.025,0.025,0.875,40)
%   tetrachoricExpectedValue(0.075,0.025,0.025,0.875,80)
%
%   tetrachoricExpectedValue(0.03,0.07,0.07,0.83,100)
%
%   References
%   ----------
%   Brown, M.B. & Benedetti, J.K.
%   On the mean and variance of the tetrachoric correlation coefficient.
%   Psychometrika (1977) 42: 347.
%   https://doi.org/10.1007/BF02293655
%
%   Author: Kristian Loewe

assert(abs(p11 + p12 + p21 + p22 - 1) < 10*eps);

assert(isvector(N) && all(N > 0) && all(mod(N,2) == 0));

%% defaults
method = 'bonett';

if ~exist('useConstraints','var') || isempty(useConstraints)
  useConstraints = 0;
end

%% if multiple N values were specified, call this function recursively
if ~isscalar(N)
  Er = zeros(size(N));
  for iN = 1:numel(N)
    Er(iN) = tetrachoricExpectedValue(p11,p12,p21,p22,N(iN),useConstraints);
  end
  return;
end

%% compute the expected value (p. 349, equation 4)
Er = 0;

lp11 = log(p11);
lp12 = log(p12);
lp21 = log(p21);
lp22 = log(p22);

if useConstraints            % use constraints as stated in the paper
  k = compK(p11,p12,p21,p22,N,useConstraints);

  ain = N*p11;
  bin = N*p12;
  cin = N*p21;

  for a = 0:ain+bin
    b = ain+bin-a;
    c = ain+cin-a;
    d = N-(a+b+c);

    r = tetrachoric(a/N,b/N,c/N,d/N, method);
    q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
    pt = q / k;
    Er = Er + r * pt;
  end

else                         % no constraints: use all tables
  k = compK(p11,p12,p21,p22,N,useConstraints);

  if p11 == p22 && p12 == p21
    for a = 0:N
      for b = 0:N-a
        for c = 0:N-(a+b)
          d = N-(a+b+c);
          if isequal([a b c d], [d c b a])
            r = tetrachoric(a/N,b/N,c/N,d/N, method);
            q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
            pt = q / k;
            Er = Er + r*pt;
          elseif (a > d || (a == d && b > c))
            continue;
          else
            r = tetrachoric(a/N,b/N,c/N,d/N, method);
            q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
            pt = q / k;
            Er = Er + 2*r*pt;
          end
        end
      end
    end

  else
    for a = 0:N
      for b = 0:N-a
        for c = 0:N-(a+b)
          d = N-(a+b+c);
          r = tetrachoric(a/N,b/N,c/N,d/N, method);
          q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
          pt = q / k;
          Er = Er + r * pt;
        end
      end
    end
  end
end

if isa(Er, 'sym')
  Er = double(Er);
end

assert(isfinite(Er));

end


function k = compK(p11,p12,p21,p22,N,useConstraints)
%COMPK Compute the coefficient k (p. 349, equation 3).
%
%   k = compK(p11,p12,p21,p22,N) computes the coefficient k based on the
%   specified cell probabilities p11,p12,p21,p22 and the total frequency N.
%
%   Author: Kristian Loewe

k = 0;

lp11 = log(p11);
lp12 = log(p12);
lp21 = log(p21);
lp22 = log(p22);

if useConstraints            % use constraints as stated in the paper

  ain = N*p11;
  bin = N*p12;
  cin = N*p21;

  for a = 0:ain+bin
    b = ain+bin-a;
    c = ain+cin-a;
    d = N-(a+b+c);

    q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
    k = k + q;
  end

else                         % no constraints: use all tables

  if p11 == p22 && p12 == p21
    for a = 0:N
      for b = 0:N-a
        for c = 0:N-(a+b)
          d = N-(a+b+c);
          if isequal([a b c d],[d c b a])
            q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
            k = k + q;
          elseif (a > d || (a == d && b > c))
            continue;
          else
            q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
            k = k + 2*q;
          end
        end
      end
    end

  else
    for a = 0:N
      for b = 0:N-a
        for c = 0:N-(a+b)
          d = N-(a+b+c);
          q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d);
          k = k + q;
        end
      end
    end
  end
end

end


function q = computeQuot(lp11,lp12,lp21,lp22,a,b,c,d)

lppow = a*lp11 + b*lp12 + c*lp21 + d*lp22;

% lpfac = log(prod(1:a)) + log(prod(1:b)) + log(prod(1:c)) + log(prod(1:d))
lpfaca = 0; for i = 1:a; lpfaca = lpfaca + log(i); end
lpfacb = 0; for i = 1:b; lpfacb = lpfacb + log(i); end
lpfacc = 0; for i = 1:c; lpfacc = lpfacc + log(i); end
lpfacd = 0; for i = 1:d; lpfacd = lpfacd + log(i); end
lpfac = lpfaca + lpfacb + lpfacc + lpfacd;

if ~isfinite(lpfac)
  if 1
    lpfac = vpa(log(prod(1:sym(a))) + log(prod(1:sym(b))) ...
      + log(prod(1:sym(c))) + log(prod(1:sym(d))));
  else
    lpfac = log(prod(1:vpa(a))) + log(prod(1:vpa(b))) ...
      + log(prod(1:vpa(c))) + log(prod(1:vpa(d)));
  end
  lpfac = double(lpfac);
end

diff = lppow - lpfac;
if diff > -700
  q = exp(diff);
else
  q = exp(vpa(diff));
end

end
