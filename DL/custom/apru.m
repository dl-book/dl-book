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

function [D, X, errs] = apru(Y, D, lam, nit, s)

% Dictionary learning with a 1-norm penalty, Algorithm 5.1
% Algorithms APrU, version OS (Sadeghi et al, SPL 2014)
% Input:
%   Y        - signal matrix
%   D        - initial dictionary
%   lam      - tradeoff parameter for 1-norm penalty
%   nit      - number of DL iterations
%   s        - sparsity level for initialization
% Output:
%   D        - learned dictionary
%   X        - final representation matrix
%   errs     - RMSE values for all iterations

% BD 21.08.2017

if nargin < 5
  s = 2;
end

[m,N] = size(Y);
[m,n] = size(D);

X = omp(Y,D,s);  % initialize representations with OMP

% main DL loop
errs = zeros(1,nit+5);
for k = 1 : nit
  E = Y - D*X;
  errs(k) = norm(E, 'fro')/ sqrt(m*N);
  for j = 1 : n     % optimize current atom and representation
    E = E + D(:,j)*X(j,:);  % modify error
    X(j,:) = wthresh(E'*D(:,j), 's', lam/2);
    d = E*X(j,:)';
    d = d / norm(d);
    D(:,j) = d;
    E = E - d*X(j,:);
  end
end

% final refinement with fixed nonzero positions
for k = 1 : 5
  E = Y - D*X;
  for j = 1 : n     % optimize current atom and representation
    supp = find(X(j,:));
    if ~isempty(supp)
      x = X(j,supp);
      F = E(:,supp) + D(:,j)*x;  % modify error
      d = F*x';
      d = d / norm(d);
      D(:,j) = d;
      X(j,supp) = F'*d;
      E(:,supp) = F - d*X(j,supp);
    end
  end
  errs(nit+k) = norm(E, 'fro')/ sqrt(m*N);
end

