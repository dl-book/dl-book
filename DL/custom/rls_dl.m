% Copyright (c) 2017-2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

function [D, X, err] = rls_dl(Y, D, s, iternum, ff)

% RLS DL (Skretting & Engan, 2010)
% Input:
%   Y        - signal matrix
%   D        - initial dictionary
%   s        - desired sparsity level
%   iternum  - number of DL iterations
%   ff       - forgetting factor
% Output:
%   D        - learned dictionary
%   X        - final representation matrix
%   err      - RMSE values for all iterations

% BD 26.05.2017

% an option that should be input parameter and dictates the implementation
method = 1;   % 0 - standard implementation from Skretting & Engan 2010. Unstable
              % 1 - Cholesky factorization of A
              % 2 - Cholesky factorization of invA

% prepare
[m,n] = size(D);
[~,N] = size(Y);
ff = 1/ff;          % we need only the inverse of the forgetting factor

err = zeros(iternum,1);
switch method
  case 0
    invA = 1e5*eye(n);       % initialize with big diagonal matrix
  case 1
    R = eye(n);
    opts1.UT = true;
    opts1.TRANSA = true;
    opts2.UT = true;
end

for t = 1 : iternum
  y = Y(:,rem(t-1,N)+1);       % artificially extract current signal
  x = omp(y,D,s);   % compute current representation
  r = y - D*x;      % residual
  
  switch method
    case 0
      u = ff*invA*x;
      a = 1 / (1 + x'*u);
      invA = ff*invA - a*u*u';  % update inverse
    case 1
      u = linsolve(R, x, opts1);
      u = ff * linsolve(R, x, opts2);
      a = 1 / (1 + x'*u);
      R = cholupdate(R,x);      % update Cholesky factorization
  end
  
  D = D + a*r*u';
  D = normc(D);
  
  X = omp(Y,D,s);
  err(t) = norm(Y - D*X, 'fro')/sqrt(m*N);
end
