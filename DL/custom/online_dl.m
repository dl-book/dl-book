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

function [D, X, err] = online_dl(Y, D, s, iternum, ff)

% Online DL 
% Input:
%   Y        - signal matrix (Mairal et al, 2010)
%   D        - initial dictionary
%   s        - desired sparsity level
%   iternum  - number of DL iterations
%   ff       - forgetting factor
% Output:
%   D        - learned dictionary
%   X        - final representation matrix
%   err      - RMSE values for all iterations

% BD 26.05.2017

% prepare
[m,n] = size(D);
[~,N] = size(Y);

err = zeros(iternum,1);
A = 1*eye(n);       % initialize with diagonal matrix to avoid zero diagonal elements in first iterations
B = 1*D;

for t = 1 : iternum
  y = Y(:,rem(t-1,N)+1);       % artificially extract current signal
  x = omp(y,D,s);   % compute current representation
  A = ff*A + x*x';  % update product matrices
  B = ff*B + y*x';
  for j = 1 : n     % update atom j
    D(:,j) = D(:,j) + (B(:,j) - D * A(:,j)) / A(j,j);
    D(:,j) = D(:,j) / norm(D(:,j));
  end

  X = omp(Y,D,s);   % compute all representations just to be able to update the error
  err(t) = norm(Y - D*X, 'fro')/sqrt(m*N);

end
