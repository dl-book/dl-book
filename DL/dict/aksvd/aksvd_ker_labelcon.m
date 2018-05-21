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

function [A, W, V, X, errs] = aksvd_ker_labelcon(Y, A, H, W, Q, V, alpha, beta, s, kerfun, iternum)
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

% BD 21.03.2018

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
meig = min(eig(K));
if meig < 0
  K = K - 2*meig*eye(N); % better add a positive constant...
end

% DL iterations
errs = zeros(iternum,1);
for iter = 1 : iternum
  % sparse coding
  X = omp_ker_labelcon(K, A, H, W, Q, V, alpha, beta, s, []);
  
  % initialize classifier matrix with LS solution !!!
  %g = 0.1;
  %W = H*X' / (X*X' + g*eye(n));
  
  % dictionary update
  for j = 1:size(A,2)
    [~, data_indices, x] = find(X(j,:));
    % optimize dictionary atom
    if (isempty(data_indices))
      a = randn(N,1);
      w = randn(c,1);
      v = randn(n,1);
      norm_aw = sqrt(a'*K*a + alpha*w'*w + beta*v'*v);
      A(:,j) = a / norm_aw;     % normalize new atom
      W(:,j) = w / norm_aw;
      V(:,j) = v / norm_aw;
      continue;
    end
    smallX = X(:,data_indices);
    A(:,j) = 0;
    a = X(j,:)' - A*(smallX*x');   % new atom
    % optimize classifier column
    W(:,j) = 0;
    w = H(:,data_indices)*x' - W*(smallX*x');
    % optimize label transformer
    V(:,j) = 0;
    v = Q(:,data_indices)*x' - V*(smallX*x');
    % normalize all current columns together!
    norm_aw = sqrt(a'*K*a + alpha*w'*w + beta*v'*v);
    a = a / norm_aw;
    w = w / norm_aw;
    v = v / norm_aw;
    % optimize representation
    u = a'*K;
    x1 = u(data_indices) - (u*A)*smallX;  % contribution of the atom
    x2 = w'*H(:,data_indices) - (w'*W)*smallX;
    x3 = v'*Q(:,data_indices) - (v'*V)*smallX;
    X(j,data_indices) = (x1+alpha*x2+beta*x3) / (1+alpha+beta);
    % store updates
    A(:,j) = a;
    W(:,j) = w;
    V(:,j) = v;
  end
  
  % compute error
  errs(iter) = norm(eye(N) - A*X,'fro') / sqrt(N*n);
end

% norm atoms to 1 and modify norms of classifier matrix
for j = 1:n
  norm_a = sqrt(A(:,j)'*K*A(:,j));
  A(:,j) = A(:,j) / norm_a;
  W(:,j) = W(:,j) / norm_a;
  V(:,j) = V(:,j) / norm_a;
  X(j,:) = X(j,:) * norm_a;
end

end

%%----------------------------------------------------------------------
function [X, shared] = omp_ker_labelcon(K, A, H, W, Q, V, alpha, beta, s, shared, varargin)
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

% BD 21.03.2018

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
        X = omp_sparse(A'*K + alpha*W'*H + beta*V'*Q, A'*K*A + alpha*W'*W + beta*V'*V, s, ompparams{:});
    else
        X = omp_err(A'*K + alpha*W'*H + beta*V'*Q, sum(K + alpha*H'*H + beta*V'*V), A'*K*A + alpha*W'*W + beta*V'*V, error, ompparams{:});
    end
end
