% Copyright (c) 2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

function [A, X, errs] = aksvd_ker(Y, A, s, kerfun, iternum)
% Kernel AK-SVD (as in Algorithm 9.2)
% Input:
%   Y        - signal matrix
%   A        - initial dictionary
%   s        - sparsity level
%   kerfun   - handle of kernel function
%   iternum  - number of DL iterations
% Output:
%   A        - learned dictionary
%   X        - final representation matrix
%   errs     - RMSE values for all iterations

% BD 24.12.2017

[N,n] = size(A);
m = size(Y,1);

% compute kernel matrix
K = zeros(N,N);
for i = 1:N
  for j = 1:i
    K(i,j) = kerfun(Y(:,i), Y(:,j));
    K(j,i) = K(i,j);
  end
end

% DL iterations
errs = zeros(iternum,1);
for iter = 1 : iternum
  % sparse coding
  X = omp_ker(K, [], A, s, []);
  
  % dictionary update
  for j = 1:size(A,2)     
    [~, data_indices, x] = find(X(j,:));
    a = A(:,j);
    A(:,j) = 0;
    if (isempty(data_indices))
      % D(:,j) = new_atom(replatoms,Y,D,X,j);
      a = randn(N,1);
      A(:,j) = a / sqrt(a'*K*a);     % normalize new atom
      continue;
    end
    smallX = X(:,data_indices);
    a = X(j,:)' - A*(smallX*x');   % new atom
    a = a / sqrt(a'*K*a);     % normalize new atom
    v = a'*K;
    X(j,data_indices) = v(data_indices) - v*A*smallX;
    A(:,j) = a;
  end
  
  % compute error
  errs(iter) = norm(eye(N) - A*X,'fro') / sqrt(N*n);
end