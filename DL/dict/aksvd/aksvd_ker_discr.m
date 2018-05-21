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

function [A, W, X, errs] = aksvd_ker_discr(Y, A, H, W, alpha, s, kerfun, iternum)
% Kernel AK-SVD for discriminative DL as in Section 9.6.3
% Input:
%   Y        - signal matrix
%   A        - initial dictionary
%   H        - label matrix
%   W        - initial classifier matrix
%   alpha    - tradeoff parameter for label error
%   s        - sparsity level
%   kerfun   - handle of kernel function
%   iternum  - number of DL iterations
% Output:
%   A        - learned dictionary
%   W        - classifier matrix
%   X        - final representation matrix
%   errs     - RMSE values for all iterations

% BD 18.03.2018

[N,n] = size(A);
%m = size(Y,1);
c = size(H,1);  % number of classes

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
  X = omp_ker_discr(K, A, H, W, alpha, s, []);
  
  % dictionary update
  for j = 1:size(A,2)
    [~, data_indices, x] = find(X(j,:));
    % optimize dictionary atom
    a = A(:,j);
    A(:,j) = 0;
    if (isempty(data_indices))
      % D(:,j) = new_atom(replatoms,Y,D,X,j);
      a = randn(N,1);
      A(:,j) = a / sqrt(a'*K*a);     % normalize new atom
      w = randn(c,1);
      W(:,j) = w/norm(w);
      continue;
    end
    smallX = X(:,data_indices);
    a = X(j,:)' - A*(smallX*x');   % new atom
    a = a / sqrt(a'*K*a);     % normalize new atom
    % optimize classifier column
    w = W(:,j);
    W(:,j) = 0;
    w = H(:,data_indices)*x' - W*(smallX*x');
    w = w / norm(w);
    % optimize representation
    v = a'*K;
    x1 = v(data_indices) - (v*A)*smallX;  % contribution of the atom
    x2 = w'*H(:,data_indices) - (w'*W)*smallX;
    X(j,data_indices) = (x1+alpha*x2) / (1+alpha);
    % store updates
    A(:,j) = a;
    W(:,j) = w;
  end
  
  % compute error
  errs(iter) = norm(eye(N) - A*X,'fro') / sqrt(N*n);
end

%%----------------------------------------------------------------------
function [X, shared] = omp_ker_discr(K, A, H, W, alpha, s, shared, varargin)
% Kernel Orthogonal Matching Pursuit algorithm, only for use inside
% kernel discriminative DL algorithm (aksvd_ker_discr)
% Input:
%   K       - kernel matrix
%   A       - current dictionary
%   H       - label matrix
%   W       - current classifier matrix
%   alpha   - trade-off parameter  
%   s       - sparsity constraint
%
% Output:
%   X       - sparse representations

% BD 18.03.2018

    if isempty(varargin)
        ompparams = {'checkdict', 'off'};
        do_sparse = true;
    else
        if strcmp(varargin{1}, 'error')
            do_sparse = false;
            error = varargin{2};
            if length(varargin) > 2
                ompparams = varargin{3:end};
            else
                ompparams = {'checkdict', 'off'};
            end
        else
            do_sparse = true;
            ompparams = varargin;
        end
    end
    if do_sparse
        X = omp_sparse(A'*K + alpha*W'*H, A'*K*A + alpha*W'*W, s, ompparams{:});
    else
        X = omp_err(A'*K + alpha*W'*H, sum(K + alpha*H'*H), A'*K*A + alpha*W'*W, error, ompparams{:});
    end
